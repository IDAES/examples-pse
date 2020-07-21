Example from:

Wilson, Zachary T., and Nikolaos V. Sahinidis. "Automated learning of
chemical reaction networks." Computers & Chemical Engineering 127
(2019): 88-98. https://doi.org/10.1016/j.compchemeng.2019.05.020

*Case 1: Isothermal CSTR*

For isothermal CSTRs across a known range of feed concentrations,
:math:`C_s^l \leq C_s^0 \leq C_s^u`, s :math:`\in` F.

The simulated reaction networks id defined below, where
:math:`k_1^{true} = 1.5`, :math:`k_2^{true} = 2.1`, and
:math:`k_3^{true} = 0.9` with a residence time of :math:`\tau = 1` is
used for the reactor.

:math:`A + B \rightarrow C \quad \{{k_1^{true}}\}`

:math:`B + C \rightarrow D \quad \{{k_2^{true}}\}`

:math:`A + D \rightarrow E \quad \{{k_3^{true}}\}`

Initial concentrations are specificed for species :math:`F = {A,B}` over
the range :math:`0 \leq C_s^0 \leq 10`, :math:`s\in F`.

.. code:: ipython3

    # Imports
    
    import pyomo.environ as pyo
    from idaes.surrogate import ripe
    import numpy as np
    import random
    import isotsim
    
    np.random.seed(20)
    
    # Setup the problem
    noise = 0.1
    ns = 5  # number of species
    lb_conc = [0,0,0,0,0]
    ub_conc = [10,10,0,0,0]
    
    # initial concentrations - only 2 data points
    cdata0 = [[1,1,0,0,0],[10,10,0,0,0]]
    cdata = isotsim.sim(cdata0)
    nd = len(cdata0) # number of data points
    
    # Expected variance based off the noise in the data
    sigma = np.multiply(noise**2,np.array(cdata))

The postulated set reaction stoichiometries is defined as:

$ A + B :raw-latex:`\rightarrow `C $

$ B + C :raw-latex:`\rightarrow `D $

$ A + D :raw-latex:`\rightarrow `E $

$ A + 2B :raw-latex:`\rightarrow `D $

$ 2A + 2B :raw-latex:`\rightarrow `E $

$ A + B + C :raw-latex:`\rightarrow `E $

$ 2A + B + D :raw-latex:`\rightarrow `C + E$

$ C + D :raw-latex:`\rightarrow `E + A $

$ 3A + 3B :raw-latex:`\rightarrow `C + E $

$ 3A + 4B :raw-latex:`\rightarrow `D + E $

$ 2A + 3B :raw-latex:`\rightarrow `C + D $

$ 4A + 5B :raw-latex:`\rightarrow `C + D + E $

.. code:: ipython3

    # considered reaction stoichiometries
    #            A   B   C   D   E
    stoich = [[ -1, -1,  1,  0,  0],
              [  0, -1, -1,  1,  0],
              [ -1,  0,  0, -1,  1],
              [ -1, -2,  0,  1,  0],
              [ -2, -2,  0,  0,  1],
              [ -1, -1, -1,  0,  1],
              [ -2, -1,  1, -1,  1],
              [  1,  0, -1, -1,  1],
              [ -3, -3,  1,  0,  1],
              [ -3, -4,  0,  1,  1],
              [ -2, -3,  1,  1,  0],
              [ -4, -5,  1,  1,  1]]

We have the initial conditions and possible stoichiometries to consider,
but we still need the kinetics reaction mechanisms. Reaction mechanisms
require a stoichiometry and kinetic model. In this case, we will be
using mass action kinetics for all the stoichiometries available, which
is built into RIPE.

.. code:: ipython3

    # IRIPE internal mass action kinetics are specified
    rxn_mechs = [['all','massact']]

Now we are ready to run the RIPE model builder:

.. code:: ipython3

    results = ripe.ripemodel(cdata,
                                 stoich = stoich,
                                 mechanisms=rxn_mechs,
                                 x0=cdata0,
                                 hide_output=False,
                                 sigma=sigma,
                                 deltaterm=0,
                                 expand_output=True)


.. parsed-literal::

    /home/ksb/anaconda3/envs/examples-rel/lib/python3.7/site-packages/numpy/core/_asarray.py:83: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray
      return array(a, dtype, copy=False, order=order)


.. parsed-literal::

       ---- Calculating null values for model selection ----    
     - Null model BIC = 7676.227335597764
     - Solving RIPE model with cardinality constraint = 1 - 
     - 1-term model BIC = 217.65886843139467
      ----  RIPE selected 1 models in the optimal reaction network  ----  
     -+-  Overall R2 of selected model : 0.98250279281449
     - Stochiometry selected for reaction 1
             [-4, -5, 1, 1, 1]
       Mechanism selected for reaction 1
             massact
       Kinetic rate constant & 95% confidence interval
             0.005408483826394419 +/- 0.0004545044904325789
    


Based on the number of data points, the best model chosen is only one
reaction.

$ 4A + 5B :raw-latex:`\rightarrow `C + D + E $

So similar to how ALAMO iterates between developing a model and adding
additional points with error maximization, RIPE provides methods to do
the same.

.. code:: ipython3

    # Adaptive experimental design using error maximization sampling
    [new_points, err] = ripe.ems(results,
                                 isotsim.sim,
                                 lb_conc,
                                 ub_conc,
                                 5, #number of species
                                 x=cdata,
                                 x0=cdata0) 
    print("New Point", new_points)
    print("Error", err)
    
    # Implement EMS as described in the RIPE publication
    new_res = isotsim.sim(new_points)[0]
    print("New Result",new_res)


.. parsed-literal::

    New Point [1.616002264907821, 0.8597699700036876, 0, 0, 0]
    Error [0.40062381 0.52402287 0.25033595 0.07038126 0.070624  ]
    New Result [1.1802669442876992, 0.22957450138430827, 0.28001096725774693, 0.0906935593637621, 0.08524270008546224]


Running ripe.ems gives us additional points maximizing the error that we
can use to develop a new model until our error tolerance is achieved. A
common loop in using RIPE follows:

.. code:: ipython3

    ite = 0
    
    while any(err >  [2*noise*s for s in new_res] ):
        print('Which concentrations violate error (True=violation) : ', err > [noise*s for s in new_res])
        results = {}
        ite+=1
        
        # Data updated explicitly 
        # so RBFopt subroutines produce consistent results
        
        new_cdata0 = np.zeros([nd+ite,ns])
        new_cdata  = np.zeros([nd+ite,ns])
        new_cdata0[:-1][:] = cdata0[:][:]
        new_cdata[:-1][:] = cdata[:][:]
        new_cdata0[-1][:] = new_points[:]
        res = isotsim.sim(new_points)[0]
        for j in range(len(res)):
            new_cdata[-1][j] = res[j]
    
        #Update weight parameters
        sigma =  np.multiply(noise**2,np.array(new_cdata))
    
        # Build updated RIPE model
        results = ripe.ripemodel(new_cdata, 
                                 stoich = stoich,
                                 mechanisms=rxn_mechs,
                                 x0=new_cdata0,
                                 sigma=sigma,
                                 expand_output=True)
    
        # Another call to EMS
        [new_points, err] = ripe.ems(results,
                                     isotsim.sim,
                                     lb_conc,
                                     ub_conc,
                                     5,
                                     x=cdata,
                                     x0=cdata0)
    
        # Update results
        new_res = isotsim.sim(new_points)[0]
        cdata0 = new_cdata0
        cdata = new_cdata


.. parsed-literal::

    Which concentrations violate error (True=violation) :  [ True  True  True  True  True]
       ---- Calculating null values for model selection ----    
     - Null model BIC = 7923.4746650041425
     - Solving RIPE model with cardinality constraint = 1 - 
     - 1-term model BIC = 465.29526828823333
     - Solving RIPE model with cardinality constraint = 2 - 
     - 2-term model BIC = 247.8231901116429
      ----  RIPE selected 2 models in the optimal reaction network  ----  
     -+-  Overall R2 of selected model : 0.9780623549758147
     - Stochiometry selected for reaction 1
             [-1, -1, 1, 0, 0]
       Mechanism selected for reaction 1
             massact
       Kinetic rate constant & 95% confidence interval
             1.2276326375145985 +/- 0.1526393150802259
    
     - Stochiometry selected for reaction 2
             [-1, -1, -1, 0, 1]
       Mechanism selected for reaction 2
             massact
       Kinetic rate constant & 95% confidence interval
             0.7397644388507576 +/- 0.09994203913081406
    
    Which concentrations violate error (True=violation) :  [ True False  True  True  True]
       ---- Calculating null values for model selection ----    
     - Null model BIC = 10404.196133123156
     - Solving RIPE model with cardinality constraint = 1 - 
     - 1-term model BIC = 1377.7742091557573
     - Solving RIPE model with cardinality constraint = 2 - 
     - 2-term model BIC = 462.64405383535995
     - Solving RIPE model with cardinality constraint = 3 - 
     - 3-term model BIC = 96.46959417194722
      ----  RIPE selected 3 models in the optimal reaction network  ----  
     -+-  Overall R2 of selected model : 0.9919042874471127
     - Stochiometry selected for reaction 1
             [-1, -1, 1, 0, 0]
       Mechanism selected for reaction 1
             massact
       Kinetic rate constant & 95% confidence interval
             1.3541627569237824 +/- 0.35391832726626726
    
     - Stochiometry selected for reaction 2
             [0, -1, -1, 1, 0]
       Mechanism selected for reaction 2
             massact
       Kinetic rate constant & 95% confidence interval
             2.0167317361130337 +/- 0.6747005426635322
    
     - Stochiometry selected for reaction 3
             [-1, 0, 0, -1, 1]
       Mechanism selected for reaction 3
             massact
       Kinetic rate constant & 95% confidence interval
             0.7702814796423452 +/- 0.28583976880203665
    


The results can vary, but RIPE can identify the simulated system of:

:math:`A + B \rightarrow C \quad \{{k_1^{true}}\}`

:math:`B + C \rightarrow D \quad \{{k_2^{true}}\}`

:math:`A + D \rightarrow E \quad \{{k_3^{true}}\}`

.. code:: ipython3

    # Final call to RIPE to get concise output
    results = ripe.ripemodel(cdata,
                             stoich = stoich,
                             mechanisms=rxn_mechs,
                             x0=cdata0,
                             sigma=sigma,
                             expand_output=False)
    print(results)


.. parsed-literal::

       ---- Calculating null values for model selection ----    
     - Null model BIC = 10404.196133123156
     - Solving RIPE model with cardinality constraint = 1 - 
     - 1-term model BIC = 1377.7742091557573
     - Solving RIPE model with cardinality constraint = 2 - 
     - 2-term model BIC = 462.64405383535995
     - Solving RIPE model with cardinality constraint = 3 - 
     - 3-term model BIC = 96.46959417194722
      ----  RIPE selected 3 models in the optimal reaction network  ----  
     -+-  Overall R2 of selected model : 0.9919042874471127
     - Stochiometry selected for reaction 1
             [-1, -1, 1, 0, 0]
       Mechanism selected for reaction 1
             massact
       Kinetic rate constant & 95% confidence interval
             1.3541627569237824 +/- 0.35391832726626726
    
     - Stochiometry selected for reaction 2
             [0, -1, -1, 1, 0]
       Mechanism selected for reaction 2
             massact
       Kinetic rate constant & 95% confidence interval
             2.0167317361130337 +/- 0.6747005426635322
    
     - Stochiometry selected for reaction 3
             [-1, 0, 0, -1, 1]
       Mechanism selected for reaction 3
             massact
       Kinetic rate constant & 95% confidence interval
             0.7702814796423452 +/- 0.28583976880203665
    
    {'k': [1.3541627569237824, 2.0167317361130337, 0.7702814796423452], 'conf_inv': [0.35391832726626726, 0.6747005426635322, 0.28583976880203665], 'stoichiometry': [[-1, -1, 1, 0, 0], [0, -1, -1, 1, 0], [-1, 0, 0, -1, 1]], 'mechanisms': ['massact', 'massact', 'massact']}


