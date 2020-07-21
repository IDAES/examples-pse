Example from:

Wilson, Zachary T., and Nikolaos V. Sahinidis. "Automated learning of
chemical reaction networks." Computers & Chemical Engineering 127
(2019): 88-98. https://doi.org/10.1016/j.compchemeng.2019.05.020

Case 2: Dynamic Chemical Looping Combusion Reactor

This is an example of a CLC reactor. The kinetic reaction rates
encapsulate solid-gas reactions. The kinetic rate laws for this example
are semi-physical or empirical to provide insights on the underlying
physical mechanisms.

The rate laws are often expressed in terms similar to

$ :raw-latex:`\frac{dX}{dT}` = kA(X)g(F) $

where $ A(X) $ is a mechanism-dependent activity term, and the function
$ g(F) $ is a parametric function of processss conditions, typically in
this case the partial pressure of methane used as a fuel.

.. code:: ipython3

    # Imports and data
    
    from idaes.surrogate import ripe
    import numpy as np
    from idaes.surrogate.ripe import mechs as mechs
    
    
    np.random.seed(20)
    
    # Import data from csv
    data = np.genfromtxt('clc.csv', delimiter=',')
    t = data[:,0]
    xdata = data[:,1]
    
    # Stoichiometry 
    # One species
    stoich = [1]


We are going to use empirical pre-defined functions from RIPE, defined
in the idaes.surrogate.ripe.mechs. The mechanisms depend on only one
species for these pre-defined rate forms.

+-------+-------+
| Rate  | $     |
| Equat | A(x)  |
| ion   | $     |
+=======+=======+
| Rando | $ 1-x |
| m     | $     |
| Nucle |       |
| ation |       |
+-------+-------+
| Power | $     |
| law   | nx^{( |
| $n =  | n-1/n |
| 2/3,  | )}    |
| 1.5,  | $     |
| 2, 3, |       |
| 4 $   |       |
+-------+-------+
| Avram | $     |
| i-Ero | n(1-x |
| feev  | )(-lo |
| $ n = | g(1-x |
| 0.5,  | ))^{( |
| 1.5,  | n-1/n |
| 2, 3, | )}    |
| 4$    | $     |
+-------+-------+
| Prout | $     |
| Tompk | x(x-1 |
| ins   | )     |
|       | $     |
+-------+-------+
| Jande | $     |
| r     | 3(1-x |
|       | )^{1/ |
|       | 3}    |
|       | (1/(1 |
|       | +x)^{ |
|       | ((-1/ |
|       | 3)-1) |
|       | })    |
|       | $     |
+-------+-------+
| Antij | $     |
| ander | 3/2(1 |
|       | -x)\  |
|       | :sup: |
|       | `{(2/ |
|       | 3)}(1 |
|       | /(1+x |
|       | )`\ { |
|       | ((-1/ |
|       | 3)-1) |
|       | })    |
|       | $     |
+-------+-------+
| Valen | $     |
| si    | 1/(-l |
|       | og(1- |
|       | x))   |
|       | $     |
+-------+-------+
| Parab | $     |
| olic  | 1/2x  |
|       | $     |
+-------+-------+
| Ginst | $     |
| ling- | (3/2) |
| Broun | (1-x) |
| tstei | :sup: |
| n     | `{(4/ |
| diffu | 3)}/( |
| sion- | (1-x) |
| 3d    | `\ {( |
|       | -1/3) |
|       | }-1)  |
|       | $     |
+-------+-------+
| Zhura | $     |
| lev-L | (3/2) |
| eshok | /((1- |
| in-Te | x)^{( |
| mpelm | -1/3) |
| an    | }-1)  |
|       | $     |
+-------+-------+
| Grain | $     |
| model | (1-x) |
|       | ^{(2/ |
|       | 3)}   |
|       | $     |
+-------+-------+

.. code:: ipython3

    # User pre-defined clc rate forms found in RIPE
    # mechs = ripe.clcforms
    clc_mechs = [mechs.randomnuc, 
                 mechs.powerlawp5, 
                 mechs.powerlaw2, 
                 mechs.powerlaw3, 
                 mechs.powerlaw4, 
                 mechs.avrami2, 
                 mechs.avrami3, 
                 mechs.avrami4, 
                 mechs.avrami5, 
                 mechs.ptompkins, 
                 mechs.jander, 
                 mechs.antijander, 
                 mechs.valensi,
                 mechs.parabolic, 
                 mechs.gb3d, 
                 mechs.zlt, 
                 mechs.grain]

All that is left is to run the ripe modeler:

.. code:: ipython3

    # Identify optimal kinetic mechanism
    results = ripe.ripemodel(xdata,
                             stoichiometry=stoich, 
                             mechanisms=clc_mechs, 
                             time=t)



.. parsed-literal::

    Calcuating residual variance with data provided with cc = 1
       ---- Calculating null values for model selection ----    
     - Null model BIC = 1852.1442125511533
     - Solving RIPE model with cardinality constraint = 1 - 
     - 1-term model BIC = 35.49650756146648
      ----  RIPE selected 1 models in the optimal reaction network  ----  
     -+-  Overall R2 of selected model : 0.18192441981887375
     - Stochiometry selected for reaction 1
             [1]
       Mechanism selected for reaction 1
             <function antijander at 0x7f57b5e87290>
       Kinetic rate constant & 95% confidence interval
             0.33776767685854847 +/- 44.373722881991505
    


