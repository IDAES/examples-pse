Module 4b: Parameter Estimation Using the NRTL State Block
----------------------------------------------------------

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

    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:29:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.016
    Total CPU secs in NLP function evaluations           =      0.014
    
    EXIT: Optimal Solution Found.


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

    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.0840425e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.8562810e-02 1.31e+03 4.38e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  6.8777957e-03 5.53e+01 2.68e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  7.1281927e-03 8.38e+01 6.43e-04  -1.7 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  7.0235071e-03 1.37e+01 1.20e-04  -1.7 4.44e-02    -  1.00e+00 1.00e+00h  1
       5  6.9830128e-03 1.13e+00 1.03e-04  -2.5 1.26e-02    -  1.00e+00 1.00e+00h  1
       6  6.7234781e-03 1.05e+02 9.56e-04  -3.8 1.18e-01    -  1.00e+00 1.00e+00h  1
       7  6.6321351e-03 5.33e-03 4.64e-05  -3.8 4.60e-03  -2.0 1.00e+00 1.00e+00h  1
       8  6.6009662e-03 4.35e-01 2.71e-05  -5.7 7.92e-03  -2.5 1.00e+00 1.00e+00h  1
       9  6.4936271e-03 5.19e+00 6.25e-05  -5.7 2.74e-02  -3.0 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.9139941e-03 1.37e+02 1.53e-03  -5.7 1.42e-01  -3.4 1.00e+00 1.00e+00h  1
      11  4.9903989e-03 7.73e+02 1.24e-02  -5.7 1.31e+00  -3.9 1.00e+00 2.50e-01h  3
      12  4.9543178e-03 1.46e+02 1.28e-03  -5.7 1.47e-01    -  1.00e+00 1.00e+00h  1
      13  4.9102927e-03 2.34e+00 3.61e-05  -5.7 1.76e-02    -  1.00e+00 1.00e+00h  1
      14  4.9095797e-03 1.27e-03 6.69e-08  -5.7 5.30e-04    -  1.00e+00 1.00e+00h  1
      15  4.9095792e-03 2.27e-05 9.72e-10  -8.6 6.27e-05    -  1.00e+00 1.00e+00h  1
      16  4.9095792e-03 2.18e-11 2.22e-14  -8.6 1.12e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 16
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9095792454458033e-03    4.9095792454458033e-03
    Dual infeasibility......:   2.2154858792359846e-14    2.2154858792359846e-14
    Constraint violation....:   5.6843418860808015e-14    2.1827872842550278e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 20
    Number of objective gradient evaluations             = 17
    Number of equality constraint evaluations            = 20
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 17
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 16
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.017
    Total CPU secs in NLP function evaluations           =      0.015
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.4165564e+01 2.91e+00 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.5264351e-02 1.21e+03 4.07e-01  -1.0 1.27e+04    -  1.00e+00 1.00e+00f  1
       2  6.0378770e-03 2.48e+02 5.72e-02  -1.7 4.05e+02    -  1.00e+00 1.00e+00h  1
       3  5.7932426e-03 2.67e+01 3.67e-04  -1.7 4.87e-01    -  1.00e+00 1.00e+00h  1
       4  5.7278926e-03 2.16e-01 5.53e-05  -2.5 6.23e-03    -  1.00e+00 1.00e+00h  1
       5  5.5451396e-03 8.13e+01 7.72e-04  -3.8 1.05e-01    -  1.00e+00 1.00e+00h  1
       6  5.4780967e-03 3.06e-03 3.62e-05  -3.8 3.59e-03  -2.0 1.00e+00 1.00e+00h  1
       7  5.4592123e-03 2.63e-01 2.11e-05  -5.7 6.17e-03  -2.5 1.00e+00 1.00e+00h  1
       8  5.3966271e-03 2.91e+00 3.76e-05  -5.7 2.06e-02  -3.0 1.00e+00 1.00e+00h  1
       9  5.1171143e-03 5.34e+01 6.12e-04  -5.7 8.85e-02  -3.4 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.3675319e-03 3.29e+03 4.77e-02  -5.7 8.03e-01  -3.9 1.00e+00 1.00e+00H  1
      11  4.1623391e-03 1.86e+03 2.54e-02  -5.7 3.91e-01    -  1.00e+00 5.00e-01h  2
      12  4.5845299e-03 1.11e+03 1.38e-02  -5.7 3.41e-01    -  1.00e+00 5.00e-01h  2
      13  4.9056458e-03 6.93e+02 8.00e-03  -5.7 2.86e-01    -  1.00e+00 5.00e-01h  2
      14  5.0852032e-03 5.36e+02 6.20e-03  -5.7 1.93e-01    -  1.00e+00 2.50e-01h  3
      15  5.1047539e-03 5.19e+02 6.00e-03  -5.7 1.63e-01    -  1.00e+00 3.12e-02h  6
      16  5.8813264e-03 1.84e+02 2.20e-03  -5.7 1.60e-01    -  1.00e+00 1.00e+00h  1
      17  5.6861869e-03 2.19e+02 1.53e-03  -5.7 1.72e-01    -  1.00e+00 1.00e+00h  1
      18  5.3975869e-03 4.57e+01 5.20e-04  -5.7 8.68e-02  -3.5 1.00e+00 1.00e+00h  1
      19  5.1548921e-03 8.20e+00 1.69e-04  -5.7 3.63e-02  -3.1 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  4.6402249e-03 1.47e+02 2.08e-03  -5.7 1.49e-01  -3.5 1.00e+00 1.00e+00h  1
      21  4.5161499e-03 1.14e+02 3.41e-03  -5.7 1.41e-01    -  1.00e+00 1.00e+00h  1
      22  4.3957175e-03 2.73e+01 5.57e-04  -5.7 6.67e-02    -  1.00e+00 1.00e+00h  1
      23  4.3886317e-03 1.29e+00 2.01e-05  -5.7 1.42e-02    -  1.00e+00 1.00e+00h  1
      24  4.3884491e-03 1.80e-03 2.42e-08  -5.7 5.20e-04    -  1.00e+00 1.00e+00h  1
      25  4.3884486e-03 3.81e-05 1.23e-09  -8.6 8.07e-05    -  1.00e+00 1.00e+00h  1
      26  4.3884486e-03 5.82e-11 2.18e-14  -8.6 1.99e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 26
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.3884486278640850e-03    4.3884486278640850e-03
    Dual infeasibility......:   2.1811060087545608e-14    2.1811060087545608e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 38
    Number of objective gradient evaluations             = 27
    Number of equality constraint evaluations            = 38
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 27
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 26
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.029
    Total CPU secs in NLP function evaluations           =      0.026
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.9510090e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.2125576e-02 1.42e+03 4.72e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.8020443e-03 7.16e+00 1.79e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.6981882e-03 1.00e+01 4.33e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  6.6662112e-03 5.68e-02 3.18e-05  -3.8 3.17e-03  -2.0 1.00e+00 1.00e+00h  1
       5  6.6220904e-03 6.43e-01 3.34e-05  -5.7 9.71e-03  -2.5 1.00e+00 1.00e+00h  1
       6  6.4622184e-03 7.93e+00 9.66e-05  -5.7 3.42e-02  -3.0 1.00e+00 1.00e+00h  1
       7  5.6129476e-03 2.26e+02 2.67e-03  -5.7 1.84e-01  -3.4 1.00e+00 1.00e+00h  1
       8  5.2786948e-03 4.84e+02 1.25e-02  -5.7 5.01e-01    -  1.00e+00 5.00e-01h  2
       9  4.9667610e-03 1.04e+02 2.05e-03  -5.7 1.28e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.9312419e-03 6.05e+00 7.46e-05  -5.7 2.99e-02    -  1.00e+00 1.00e+00h  1
      11  4.9300462e-03 1.22e-02 9.93e-08  -5.7 1.30e-03    -  1.00e+00 1.00e+00h  1
      12  4.9300439e-03 1.65e-08 8.43e-14  -5.7 1.47e-06    -  1.00e+00 1.00e+00h  1
      13  4.9300437e-03 2.20e-05 1.04e-09  -8.6 6.18e-05    -  1.00e+00 1.00e+00h  1
      14  4.9300437e-03 5.82e-11 2.09e-14  -8.6 1.08e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9300436762568598e-03    4.9300436762568598e-03
    Dual infeasibility......:   2.0897124089924014e-14    2.0897124089924014e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 17
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 17
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.013
    Total CPU secs in NLP function evaluations           =      0.012
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.9044622e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.3440372e-02 1.31e+03 4.28e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  6.7701054e-03 1.94e+01 1.56e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  6.7632185e-03 5.13e-01 8.29e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  4.9984577e-03 6.53e+03 4.99e-02  -3.8 9.56e-01    -  8.83e-01 1.00e+00h  1
       5  5.8970606e-03 3.69e+03 2.86e-02  -3.8 5.24e-01    -  1.00e+00 5.00e-01h  2
       6  6.2284432e-03 1.02e+03 5.66e-03  -3.8 3.79e-01    -  1.00e+00 1.00e+00h  1
       7  5.4575879e-03 8.42e+02 1.23e-02  -3.8 5.04e-01  -4.0 1.00e+00 5.00e-01h  2
       8  4.9500121e-03 1.78e+02 2.42e-03  -3.8 1.63e-01    -  1.00e+00 1.00e+00h  1
       9  4.8925167e-03 1.05e+01 9.72e-05  -3.8 3.84e-02    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8905079e-03 1.91e-02 9.75e-08  -3.8 1.56e-03    -  1.00e+00 1.00e+00h  1
      11  4.8891793e-03 1.59e-01 6.94e-06  -5.7 5.25e-03    -  1.00e+00 1.00e+00h  1
      12  4.8891580e-03 3.87e-05 5.52e-10  -5.7 7.76e-05    -  1.00e+00 1.00e+00h  1
      13  4.8891578e-03 2.30e-05 1.03e-09  -8.6 6.31e-05    -  1.00e+00 1.00e+00h  1
      14  4.8891578e-03 8.73e-11 2.19e-14  -8.6 1.13e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8891577511399383e-03    4.8891577511399383e-03
    Dual infeasibility......:   2.1924853700658986e-14    2.1924853700658986e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 18
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 18
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.011
    Total CPU secs in NLP function evaluations           =      0.010
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  7.8150843e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.2556389e-01 1.40e+03 5.22e-01  -1.0 1.37e+04    -  9.89e-01 1.00e+00f  1
       2  7.8835890e-03 1.90e+02 5.40e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  7.7016614e-03 7.33e+00 1.73e-04  -1.7 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  7.6321630e-03 1.64e+00 1.64e-04  -3.8 1.51e-02    -  1.00e+00 1.00e+00h  1
       5  7.0230261e-03 5.61e+02 4.64e-03  -3.8 2.74e-01  -4.0 1.00e+00 1.00e+00h  1
       6  5.3677119e-03 1.37e+03 1.93e-02  -3.8 1.57e+00  -3.6 1.00e+00 2.50e-01h  3
       7  5.5193320e-03 2.88e+02 2.16e-03  -3.8 2.06e-01    -  1.00e+00 1.00e+00h  1
       8  5.3458429e-03 2.26e+00 1.40e-04  -3.8 2.00e-02    -  1.00e+00 1.00e+00h  1
       9  5.3436454e-03 4.09e-03 8.90e-08  -3.8 9.00e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.3423609e-03 1.44e-01 5.59e-06  -5.7 5.00e-03    -  1.00e+00 1.00e+00h  1
      11  5.3423384e-03 3.09e-05 4.83e-10  -5.7 6.96e-05    -  1.00e+00 1.00e+00h  1
      12  5.3423382e-03 2.08e-05 8.29e-10  -8.6 6.03e-05    -  1.00e+00 1.00e+00h  1
      13  5.3423382e-03 4.37e-11 2.08e-14  -8.6 1.02e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.3423382226729888e-03    5.3423382226729888e-03
    Dual infeasibility......:   2.0783207275703022e-14    2.0783207275703022e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 17
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 17
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.010
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.4978996e+01 2.91e+00 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.5701986e-02 1.22e+03 3.86e-01  -1.0 1.27e+04    -  1.00e+00 1.00e+00f  1
       2  6.4515492e-03 3.34e+01 2.01e-02  -1.7 4.05e+02    -  1.00e+00 1.00e+00h  1
       3  6.4119874e-03 2.81e-01 6.38e-05  -2.5 4.87e-01    -  1.00e+00 1.00e+00h  1
       4  5.4928106e-03 1.10e+03 9.14e-03  -3.8 3.90e-01    -  9.78e-01 1.00e+00h  1
       5  4.7245659e-03 5.81e+00 1.59e-03  -3.8 4.53e-02    -  1.00e+00 1.00e+00h  1
       6  4.7093275e-03 3.12e-01 1.71e-05  -3.8 7.17e-03    -  1.00e+00 1.00e+00h  1
       7  4.7080691e-03 1.56e-01 6.33e-06  -5.7 5.19e-03    -  1.00e+00 1.00e+00h  1
       8  4.7080478e-03 3.70e-05 4.73e-10  -5.7 7.57e-05    -  1.00e+00 1.00e+00h  1
       9  4.7080476e-03 2.66e-05 1.08e-09  -8.6 6.77e-05    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7080476e-03 7.28e-11 2.20e-14  -8.6 1.33e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 10
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7080475529480475e-03    4.7080475529480475e-03
    Dual infeasibility......:   2.2033769380896482e-14    2.2033769380896482e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.009
    Total CPU secs in NLP function evaluations           =      0.008
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.4651711e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.3177302e-02 1.41e+03 4.82e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  7.2326237e-03 1.31e+01 1.73e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  7.2319314e-03 9.13e-01 1.07e-04  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.3076526e-03 1.04e+04 8.02e-02  -3.8 1.62e+01    -  2.43e-01 7.49e-02h  3
       5  6.0511278e-03 7.95e+03 6.19e-02  -3.8 7.10e-01    -  8.69e-01 2.50e-01h  3
       6  6.5259828e-03 6.10e+03 4.77e-02  -3.8 7.87e-01    -  1.00e+00 2.50e-01h  3
       7  6.6411904e-03 5.36e+03 4.19e-02  -3.8 5.28e-01    -  1.00e+00 1.25e-01h  4
       8  6.8219104e-03 4.12e+03 3.22e-02  -3.8 4.91e-01    -  1.00e+00 2.50e-01h  3
       9  6.8982769e-03 3.17e+03 2.46e-02  -3.8 4.24e-01    -  1.00e+00 2.50e-01h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  7.1171600e-03 9.72e+02 5.63e-03  -3.8 3.65e-01    -  1.00e+00 1.00e+00h  1
      11  5.6698956e-03 1.84e+00 1.04e-03  -3.8 8.35e-02  -2.0 1.00e+00 1.00e+00h  1
      12  5.3620427e-03 3.22e+02 8.27e-03  -3.8 2.29e-01    -  1.00e+00 1.00e+00h  1
      13  5.1557773e-03 6.15e+01 1.31e-03  -3.8 9.86e-02    -  1.00e+00 1.00e+00h  1
      14  5.1407341e-03 2.83e+00 4.36e-05  -3.8 2.07e-02    -  1.00e+00 1.00e+00h  1
      15  5.1403646e-03 4.24e-03 4.98e-08  -3.8 7.86e-04    -  1.00e+00 1.00e+00h  1
      16  5.1391094e-03 1.37e-01 6.35e-06  -5.7 4.87e-03    -  1.00e+00 1.00e+00h  1
      17  5.1390902e-03 2.76e-05 4.25e-10  -5.7 6.57e-05    -  1.00e+00 1.00e+00h  1
      18  5.1390900e-03 1.99e-05 9.42e-10  -8.6 5.88e-05    -  1.00e+00 1.00e+00h  1
      19  5.1390900e-03 1.75e-10 2.19e-14  -8.6 9.64e-09    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 19
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.1390900353160832e-03    5.1390900353160832e-03
    Dual infeasibility......:   2.1895785635083877e-14    2.1895785635083877e-14
    Constraint violation....:   6.4580138631362250e-13    1.7462298274040222e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 34
    Number of objective gradient evaluations             = 20
    Number of equality constraint evaluations            = 34
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 20
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 19
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.019
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.6106959e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.1134730e-02 1.30e+03 4.57e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  5.6814867e-03 5.03e+01 2.57e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  5.6783824e-03 2.16e-02 4.02e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  5.3650913e-03 1.95e+02 1.84e-03  -3.8 1.64e-01    -  1.00e+00 1.00e+00h  1
       5  5.1380210e-03 6.77e-02 7.09e-05  -3.8 6.99e-03  -2.0 1.00e+00 1.00e+00h  1
       6  5.0960925e-03 5.74e-01 3.28e-05  -5.7 9.38e-03  -2.5 1.00e+00 1.00e+00h  1
       7  4.9684299e-03 5.94e+00 1.01e-04  -5.7 3.03e-02  -3.0 1.00e+00 1.00e+00h  1
       8  4.6031822e-03 6.24e+01 1.06e-03  -5.7 9.89e-02  -3.4 1.00e+00 1.00e+00h  1
       9  4.4853812e-03 9.39e+01 2.45e-03  -5.7 1.26e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4202306e-03 1.73e+01 3.27e-04  -5.7 5.28e-02    -  1.00e+00 1.00e+00h  1
      11  4.4166175e-03 4.76e-01 7.32e-06  -5.7 8.57e-03    -  1.00e+00 1.00e+00h  1
      12  4.4165512e-03 2.23e-04 2.96e-09  -5.7 1.82e-04    -  1.00e+00 1.00e+00h  1
      13  4.4165509e-03 3.89e-05 1.27e-09  -8.6 8.15e-05    -  1.00e+00 1.00e+00h  1
      14  4.4165509e-03 5.82e-11 2.18e-14  -8.6 2.03e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4165509401182914e-03    4.4165509401182914e-03
    Dual infeasibility......:   2.1826336974431759e-14    2.1826336974431759e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.015
    Total CPU secs in NLP function evaluations           =      0.012
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.2735527e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.2548388e-02 1.41e+03 4.98e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.3264561e-03 1.64e+01 1.76e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.3755192e-03 2.83e-02 1.80e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.8929267e-03 4.05e+02 3.94e-03  -3.8 2.38e-01    -  1.00e+00 1.00e+00h  1
       5  4.2942299e-03 4.50e+02 7.28e-03  -3.8 4.11e-01  -4.0 1.00e+00 5.00e-01h  2
       6  4.2875721e-03 7.41e+01 6.78e-04  -3.8 1.04e-01    -  1.00e+00 1.00e+00h  1
       7  4.2725280e-03 1.02e+00 8.87e-06  -3.8 1.13e-02    -  1.00e+00 1.00e+00h  1
       8  4.2722578e-03 6.07e-05 4.52e-09  -3.8 1.29e-04    -  1.00e+00 1.00e+00h  1
       9  4.2706040e-03 3.00e-01 9.81e-06  -5.7 7.12e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.2705741e-03 1.54e-04 2.04e-09  -5.7 1.54e-04    -  1.00e+00 1.00e+00h  1
      11  4.2705738e-03 4.20e-05 1.43e-09  -8.6 8.45e-05    -  1.00e+00 1.00e+00h  1
      12  4.2705738e-03 2.91e-11 2.17e-14  -8.6 2.20e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.2705738341092475e-03    4.2705738341092475e-03
    Dual infeasibility......:   2.1650569535642392e-14    2.1650569535642392e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.011
    Total CPU secs in NLP function evaluations           =      0.008
    
    EXIT: Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-05-07 21:30:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  7.0875113e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.1471761e-01 1.40e+03 5.16e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.9933860e-03 8.23e+01 3.44e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  7.1409700e-03 5.09e+01 4.02e-04  -1.7 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  7.1177587e-03 1.94e+00 1.49e-04  -2.5 1.73e-02    -  1.00e+00 1.00e+00h  1
       5  6.9488418e-03 5.00e+01 5.40e-04  -3.8 8.11e-02    -  1.00e+00 1.00e+00h  1
       6  6.0289128e-03 1.73e+03 1.34e-02  -3.8 4.85e-01  -4.0 1.00e+00 1.00e+00h  1
       7  5.4883638e-03 2.40e+02 3.18e-03  -3.8 1.79e-01    -  1.00e+00 1.00e+00h  1
       8  4.9833716e-03 3.85e+02 8.85e-04  -3.8 2.61e-01    -  1.00e+00 1.00e+00H  1
       9  5.0713750e-03 1.35e+02 7.42e-04  -3.8 1.36e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.0520063e-03 1.63e+00 1.74e-05  -3.8 1.49e-02    -  1.00e+00 1.00e+00h  1
      11  5.0515050e-03 1.37e-05 2.34e-09  -3.8 9.95e-05    -  1.00e+00 1.00e+00h  1
      12  5.0501264e-03 1.75e-01 6.52e-06  -5.7 5.51e-03    -  1.00e+00 1.00e+00h  1
      13  5.0501023e-03 4.77e-05 7.10e-10  -5.7 8.63e-05    -  1.00e+00 1.00e+00h  1
      14  5.0501021e-03 2.52e-05 9.62e-10  -8.6 6.61e-05    -  1.00e+00 1.00e+00h  1
      15  5.0501021e-03 5.82e-11 2.18e-14  -8.6 1.26e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.0501021012331858e-03    5.0501021012331858e-03
    Dual infeasibility......:   2.1815511555876159e-14    2.1815511555876159e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 17
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 17
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.011
    
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
          <td>0.490504</td>
          <td>-0.418699</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.458572</td>
          <td>-0.394421</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.494748</td>
          <td>-0.421790</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.491900</td>
          <td>-0.419662</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.505956</td>
          <td>-0.430691</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.482191</td>
          <td>-0.412305</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.503123</td>
          <td>-0.428203</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.458817</td>
          <td>-0.394662</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.451992</td>
          <td>-0.389319</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.491462</td>
          <td>-0.419709</td>
        </tr>
      </tbody>
    </table>
    </div>

