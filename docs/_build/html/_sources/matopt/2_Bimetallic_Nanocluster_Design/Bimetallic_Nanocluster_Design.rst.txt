Bimetallic Nanocluster Cohesive Energy Minimization - via Labeling
==================================================================

This notebook serves as an example application of the MatOpt framework.
We consider an example optimization problem of identifying the global
energy minimum bimetallic nanocluster configuration.

This is a continuation of the example given in
***Monometallic\_Nanocluster\_Design.ipynb***. In this example, we will
show how a very similar model can be used to optimize a bimetallic
cluster by "labelling" the sites of a pre-defined monometallic cluster.

The model for cohesive energy is based on:

Yan, Z., Taylor, M. G., Mascareno, A., & Mpourmpakis, G. (2018). Size-,
Shape-, and Composition-Dependent Model for Metal Nanoparticle Stability
Prediction. *Nano Letters*, 18(4), 2696-2704.

Importing Packages
------------------

We start by importing several standard Python modules for convienience.

.. code:: ipython3

    import numpy as np
    from math import sqrt

Then, we import the MatOpt package in its entirety.

.. code:: ipython3

    from idaes.apps.matopt import *

Setting Up a Material System
----------------------------

We first identify the optimal metal-independent nanocluster shape, using
the code that wsas demonstrated in
**Monometallic\_Nanocluster\_Design.ipynb**.

.. code:: ipython3

    Lat = FCCLattice(IAD=1.0)
    Canv = Canvas()
    Canv.addLocation(np.array([0,0,0],dtype=float))
    Canv.addShells(1,Lat.getNeighbors)
    Atoms = [Atom('Cu')]
    N = 6
    m = MatOptModel(Canv,Atoms)
    Vals = [sqrt(CN) for CN in range(0,13)]
    BPs = [CN for CN in range(0,13)]
    m.addSitesDescriptor('CNRi',bounds=(0,sqrt(12)),integer=False,
                         rules=PiecewiseLinear(values=Vals,
                                               breakpoints=BPs,
                                               input_desc=m.Ci))
    m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumSites(desc=m.CNRi,
                                                        coefs=(1/(N*sqrt(12))))))
    m.addGlobalDescriptor('Size',bounds=(N,N),
                          rules=EqualTo(SumSites(desc=m.Yi)))

.. code:: ipython3

    # This step requires the CPLEX solver
    try:
        D = m.maximize(m.Ecoh,tilim=100)
    except Exception as err:
        D = None  # rest of this notebook won't do much


.. parsed-literal::

    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer Community Edition 12.9.0.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2019.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmp7_frnrqo.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 100
    CPLEX> Problem '/tmp/tmp_h5lmt8k.pyomo.lp' read.
    Read time = 0.00 sec. (0.06 ticks)
    CPLEX> Problem name         : /tmp/tmp_h5lmt8k.pyomo.lp
    Objective sense      : Maximize
    Variables            :     426  [Nneg: 1,  Fix: 1,  Box: 13,  Free: 157,
                                     Binary: 241,  General Integer: 13]
    Objective nonzeros   :       1
    Linear constraints   :     583  [Less: 515,  Greater: 13,  Equal: 55]
      Nonzeros           :    1866
      RHS nonzeros       :      86
    
    Variables            : Min LB: 0.000000         Max UB: 12.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.04811252       Max   : 12.00000       
      RHS nonzeros       : Min   : 1.000000         Max   : 1.000000       
    CPLEX> CPXPARAM_TimeLimit                               100
    CPXPARAM_MIP_Tolerances_AbsMIPGap                0
    CPXPARAM_MIP_Tolerances_MIPGap                   0
    Tried aggregator 2 times.
    MIP Presolve eliminated 171 rows and 147 columns.
    MIP Presolve modified 12 coefficients.
    Aggregator did 38 substitutions.
    Reduced MIP has 374 rows, 241 columns, and 1149 nonzeros.
    Reduced MIP has 156 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.01 sec. (1.49 ticks)
    Found incumbent of value 0.521849 after 0.08 sec. (2.09 ticks)
    Probing time = 0.00 sec. (1.36 ticks)
    Cover probing fixed 0 vars, tightened 1 bounds.
    Tried aggregator 1 time.
    MIP Presolve eliminated 2 rows and 1 columns.
    Reduced MIP has 372 rows, 240 columns, and 1144 nonzeros.
    Reduced MIP has 156 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (1.50 ticks)
    Probing time = 0.00 sec. (1.31 ticks)
    Clique table members: 880.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 2 threads.
    Root relaxation solution time = 0.01 sec. (3.19 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
    *     0+    0                            0.5218        2.1667           315.19%
          0     0        0.9861   105        0.5218        0.9861      278   88.96%
          0     0        0.6455    96        0.5218      Cuts: 68      386   23.69%
    *     0+    0                            0.5485        0.6455            17.68%
          0     0        0.6351    99        0.5485      Cuts: 82      430   15.80%
          0     0        0.6329    93        0.5485  ZeroHalf: 26      446   15.39%
          0     0        0.6327    99        0.5485  ZeroHalf: 36      467   15.35%
          0     0        0.6322    99        0.5485      Cuts: 20      485   15.26%
          0     0        0.6037    98        0.5485  ZeroHalf: 28      538   10.07%
          0     0        0.6027    97        0.5485  ZeroHalf: 14      566    9.88%
          0     0        0.6011    94        0.5485  ZeroHalf: 18      587    9.60%
          0     0        0.6005    98        0.5485  ZeroHalf: 12      609    9.48%
          0     0        0.5992    98        0.5485  ZeroHalf: 19      645    9.25%
          0     0        0.5986    98        0.5485   ZeroHalf: 5      657    9.14%
          0     0        0.5979    94        0.5485   ZeroHalf: 5      670    9.01%
          0     0        0.5977    94        0.5485   ZeroHalf: 5      677    8.97%
          0     0        0.5970    98        0.5485   ZeroHalf: 9      698    8.84%
          0     0        0.5963    98        0.5485  ZeroHalf: 15      720    8.72%
          0     0        0.5958    93        0.5485   ZeroHalf: 3      733    8.62%
          0     0        0.5957    98        0.5485   ZeroHalf: 8      754    8.60%
          0     0        0.5955    92        0.5485  ZeroHalf: 19      772    8.57%
    *     0+    0                            0.5500        0.5955             8.27%
          0     2        0.5955    92        0.5500        0.5955      772    8.27%
    Elapsed time = 0.80 sec. (99.57 ticks, tree = 0.02 MB, solutions = 3)
    
    Clique cuts applied:  13
    Implied bound cuts applied:  14
    Zero-half cuts applied:  16
    Lift and project cuts applied:  1
    Gomory fractional cuts applied:  3
    
    Root node processing (before b&c):
      Real time             =    0.80 sec. (99.33 ticks)
    Parallel b&c, 2 threads:
      Real time             =    0.04 sec. (5.00 ticks)
      Sync time (average)   =    0.02 sec.
      Wait time (average)   =    0.00 sec.
                              ------------
    Total (root+branch&cut) =    0.84 sec. (104.33 ticks)
    
    Solution pool: 3 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  5.5003296046e-01
    Solution time =    0.84 sec.  Iterations = 1088  Nodes = 17
    Deterministic time = 104.33 ticks  (123.84 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpsx0_sgi9.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.5500329604578591


We take the locations from the optimal monometallic problem to
initialize a ***Canvas*** object for the bimetallic case.

.. code:: ipython3

    Canv = Canvas()
    if D:  # may be None if CPLEX was not found
        for i in range(len(D)):
            if(D.Contents[i] is not None):
                Canv.addLocation(D.Canvas.Points[i])
    Canv.setNeighborsFromFunc(Lat.getNeighbors)

Additionally, we create a few data structures for holding bimetallic
material information. First, we make a list of multiple ***Atom***
objects that will be the building blocks of the model. Next, we specify
a dictionary with the bounds to impose on composition.

.. code:: ipython3

    Atoms = [Atom('Cu'),Atom('Ag')]
    CompBounds = {Atom('Cu'):(3,3),
                  Atom('Ag'):(3,3)}

Specifying an Optimization Model
--------------------------------

We start by creating a ***MatOptModel*** object that will hold the
information about the problem variables and constraints. At a minimum,
ever model requires a Canvas object to be defined. Additionally, the
list of building blocks and conformations that are present in the model
should be defined.

.. code:: ipython3

    m = MatOptModel(Canv,Atoms)

By default, several basic variables are pre-defined. See the first
example, ***Monometallic\_Nanocluster\_Design.ipynb*** for a description
of basic variables, expressions, and constraint rules.

To start, we inidcate that the choice to place an atom is fixed so that
each canvas site is required to have an atom. This simplifies the
problem significantly and results in a model that will seek to find the
optimal labeling of metals on the nanocluster.

.. code:: ipython3

    m.Yi.rules.append(FixedTo(1.0))

Next, we define a descriptor for the energy of bonds as a function of
properties at each site. Since the locations of the atoms are fixed, the
only decision is how to label each site as either Atom A or Atom B. This
allows us to simplify the model and compute coefficients that rely on
coordination number. In the block below, we implement the bimetallic
model for bond energy defined in Yan et al., 2018.

.. code:: ipython3

    GklCoefs = {(Atom('Cu'),Atom('Cu')):3.520,
                (Atom('Cu'),Atom('Ag')):2.112,
                (Atom('Ag'),Atom('Ag')):2.580,
                (Atom('Ag'),Atom('Cu')):3.612}
    BEijCoefs = {}
    for i in range(len(Canv)):
        CNi = sum(1 for _ in Canv.NeighborhoodIndexes[i] if _ is not None)
        for j in Canv.NeighborhoodIndexes[i]:
            if(j is not None):
                CNj = sum(1 for _ in Canv.NeighborhoodIndexes[j] if _ is not None)
                for k in Atoms:
                    for l in Atoms:
                        BEijCoefs[i,j,k,l] = GklCoefs[k,l]*1/sqrt(CNi) + GklCoefs[l,k]*1/sqrt(CNj)
    m.addBondsDescriptor('BEij',
                         rules=EqualTo(SumBondTypes(m.Xijkl,coefs=BEijCoefs)),
                         symmetric_bonds=True)

Next, we define the cohesive energy as a sum of contributions from all
BEij bond descriptors.

.. code:: ipython3

    m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumBonds(desc=m.BEij,
                                                        coefs=1/(N*sqrt(12)))))

Finally, we add constraints on the size and composition of the resulting
designs.

.. code:: ipython3

    m.addGlobalTypesDescriptor('Composition',bounds=CompBounds,
                               rules=EqualTo(SumSites(desc=m.Yik)))

Solving the Model
-----------------

Once the model is fully specified, we can optimize in light of a global
descriptor. In this example, we choose to maximize the cohesive energy
defined previously. Additionally, we can specify basic optimization
parameters such as the time limit and memory limit\* for the optimizer.

.. code:: ipython3

    D = None
    try:
        D = m.maximize(m.Ecoh,tilim=360,trelim=4096)
    except:
        print('MaOpt can not find usable solver (CPLEX or NEOS-CPLEX)')


.. parsed-literal::

    WARNING: DEPRECATED: SetProduct.set_tuple is deprecated.  Use
        SetProduct.subsets() to get the operator arguments.  (deprecated in TBD)
        (called from /home/ksb/anaconda3/envs/examples-rel/lib/python3.7/site-
        packages/idaes/apps/matopt/../matopt/opt/pyomo_modeling.py:284)
    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer Community Edition 12.9.0.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2019.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmpc8co_b12.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 360
    CPLEX> New value for upper limit on size of tree in megabytes: 4096
    CPLEX> Problem '/tmp/tmpfwvhlu5v.pyomo.lp' read.
    Read time = 0.04 sec. (0.01 ticks)
    CPLEX> Problem name         : /tmp/tmpfwvhlu5v.pyomo.lp
    Objective sense      : Maximize
    Variables            :      71  [Nneg: 1,  Fix: 2,  Free: 12,  Binary: 56]
    Objective nonzeros   :       1
    Linear constraints   :     159  [Less: 138,  Equal: 21]
      Nonzeros           :     414
      RHS nonzeros       :      57
    
    Variables            : Min LB: 0.000000         Max UB: 3.000000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.04811252       Max   : 4.064546       
      RHS nonzeros       : Min   : 1.000000         Max   : 1.000000       
    CPLEX> CPXPARAM_TimeLimit                               360
    CPXPARAM_MIP_Tolerances_AbsMIPGap                0
    CPXPARAM_MIP_Tolerances_MIPGap                   0
    CPXPARAM_MIP_Limits_TreeMemory                   4096
    Tried aggregator 2 times.
    MIP Presolve eliminated 20 rows and 15 columns.
    Aggregator did 6 substitutions.
    Reduced MIP has 133 rows, 50 columns, and 314 nonzeros.
    Reduced MIP has 50 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (0.44 ticks)
    Found incumbent of value 1.546227 after 0.00 sec. (0.80 ticks)
    Probing time = 0.00 sec. (0.23 ticks)
    Tried aggregator 1 time.
    Reduced MIP has 133 rows, 50 columns, and 314 nonzeros.
    Reduced MIP has 50 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (0.30 ticks)
    Probing time = 0.00 sec. (0.23 ticks)
    Clique table members: 258.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 2 threads.
    Root relaxation solution time = 0.00 sec. (0.31 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
    *     0+    0                            1.5462        6.5036           320.61%
    *     0+    0                            1.6517        6.5036           293.74%
          0     0        3.2518    50        1.6517        3.2518       52   96.87%
          0     0        2.0902    40        1.6517     Cuts: 107       92   26.54%
    *     0+    0                            1.6556        2.0902            26.25%
          0     0        1.8431    38        1.6556      Cuts: 21      114   11.32%
          0     0        1.8217    16        1.6556      Cuts: 16      116   10.03%
    *     0+    0                            1.6754        1.8217             8.73%
          0     0        1.8116    37        1.6754   ZeroHalf: 2      120    8.13%
          0     0        1.7448    36        1.6754   ZeroHalf: 9      122    4.14%
          0     0        cutoff              1.6754                    128     --- 
    Elapsed time = 0.06 sec. (9.19 ticks, tree = 0.01 MB, solutions = 4)
    
    Clique cuts applied:  31
    Implied bound cuts applied:  8
    Zero-half cuts applied:  13
    Lift and project cuts applied:  1
    Gomory fractional cuts applied:  1
    
    Root node processing (before b&c):
      Real time             =    0.06 sec. (9.20 ticks)
    Parallel b&c, 2 threads:
      Real time             =    0.00 sec. (0.00 ticks)
      Sync time (average)   =    0.00 sec.
      Wait time (average)   =    0.00 sec.
                              ------------
    Total (root+branch&cut) =    0.06 sec. (9.20 ticks)
    
    Solution pool: 4 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  1.6754118151e+00
    Solution time =    0.10 sec.  Iterations = 128  Nodes = 0
    Deterministic time = 9.20 ticks  (93.24 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmp6ste5uaf.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 1.6754118151174977


Processing Solutions
--------------------

If a design was identified (optimal or otherwise), then a ***Design***
object is returned from the optimization method. The optimal design can
be plotted via any of the supported parsers.

.. code:: ipython3

    if D is not None:
        D.toPDB('result.pdb')
