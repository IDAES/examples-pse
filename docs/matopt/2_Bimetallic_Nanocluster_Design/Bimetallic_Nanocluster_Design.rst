Bimetallic Nanocluster Cohesive Energy Minimization - via Labeling
==================================================================

This notebook serves as an example application of the MatOpt framework.
We consider an example optimization problem of identifying the global
energy minimum bimetallic nanocluster configuration.

This is a continuation of the example given in
**Monometallic_Nanocluster_Design.ipynb**. In this example, we will show
how a very similar model can be used to optimize a bimetallic cluster by
“labelling” the sites of a pre-defined monometallic cluster.

The model for cohesive energy is based on:

Yan, Z., Taylor, M. G., Mascareno, A., & Mpourmpakis, G. (2018). Size-,
Shape-, and Composition-Dependent Model for Metal Nanoparticle Stability
Prediction. *Nano Letters*, 18(4), 2696-2704.

Importing Packages
------------------

We start by importing several standard Python modules for convienience.

.. code:: ipython3

    import os 
    from math import sqrt

Then, we import the MatOpt package in its entirety.

.. code:: ipython3

    from matopt import *


::


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-2-be5adc2fb850> in <module>
    ----> 1 from matopt import *
    

    ModuleNotFoundError: No module named 'matopt'


Setting Up a Material System
----------------------------

We first identify the optimal metal-independent nanocluster shape, using
the code that wsas demonstrated in
**Monometallic_Nanocluster_Design.ipynb**.

.. code:: ipython3

    Lat = FCCLattice(IAD=1.0)
    Canv = Canvas()
    Canv.addLocation(np.array([0,0,0],dtype=float))
    Canv.addShells(2,Lat.getNeighbors)
    Atoms = [Atom('Cu')]
    N = 20
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
    D = m.maximize(m.Ecoh,tilim=100)


.. parsed-literal::

    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 12.6.1.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2014.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmpcvroj97s.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 100
    CPLEX> Problem '/tmp/tmpwjv0wzzl.pyomo.lp' read.
    Read time = 0.00 sec. (0.18 ticks)
    CPLEX> Problem name         : /tmp/tmpwjv0wzzl.pyomo.lp
    Objective sense      : Maximize
    Variables            :    1315  [Nneg: 716,  Fix: 1,  Box: 55,  Free: 1,
                                     Binary: 487,  General Integer: 55]
    Objective nonzeros   :       1
    Linear constraints   :    1519  [Less: 1296,  Equal: 223]
      Nonzeros           :    5769
      RHS nonzeros       :     488
    SOS                  :      55  [SOS2: 55, 715 members, all continuous]
    
    Variables            : Min LB: 0.000000         Max UB: 20.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.01443376       Max   : 12.00000       
      RHS nonzeros       : Min   : 1.000000         Max   : 1.000000       
    CPLEX> Tried aggregator 2 times.
    MIP Presolve eliminated 2 rows and 189 columns.
    Aggregator did 55 substitutions.
    Reduced MIP has 1462 rows, 1071 columns, and 5043 nonzeros.
    Reduced MIP has 487 binaries, 0 generals, 55 SOSs, and 0 indicators.
    Presolve time = 0.01 sec. (5.87 ticks)
    Probing time = 0.01 sec. (5.75 ticks)
    Tried aggregator 1 time.
    MIP Presolve eliminated 648 rows and 216 columns.
    Reduced MIP has 814 rows, 855 columns, and 3531 nonzeros.
    Reduced MIP has 271 binaries, 0 generals, 55 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (4.18 ticks)
    Probing time = 0.00 sec. (1.86 ticks)
    Tried aggregator 1 time.
    Reduced MIP has 814 rows, 855 columns, and 3531 nonzeros.
    Reduced MIP has 271 binaries, 0 generals, 55 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (2.36 ticks)
    Probing time = 0.00 sec. (1.87 ticks)
    Clique table members: 432.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 8 threads.
    Root relaxation solution time = 0.02 sec. (36.82 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
          0     0        1.3207   271                      1.3207     1079         
    *     0+    0                            0.5260        1.3207           151.11%
    *     0+    0                            0.5412        1.3207           144.03%
          0     0        0.7995   326        0.5412     Cuts: 303     1668   47.73%
    *     0+    0                            0.6097        0.7995            31.13%
          0     0        0.7922   323        0.6097       Cuts: 5     1760   29.93%
          0     0        0.7909   324        0.6097       Cuts: 6     1809   29.72%
          0     0        0.7903   325        0.6097   ZeroHalf: 7     1848   29.62%
          0     0        0.7897   326        0.6097   ZeroHalf: 4     1881   29.52%
    *     0+    0                            0.6771        0.7897            16.64%
          0     0        0.7893   326        0.6771   ZeroHalf: 5     1910   16.56%
          0     0        0.7875   319        0.6771       Cuts: 9     1969   16.31%
          0     0        0.7874   326        0.6771   ZeroHalf: 2     1989   16.29%
          0     0        0.7873   326        0.6771   ZeroHalf: 4     2001   16.28%
    *     0+    0                            0.6792        0.7873            15.92%
          0     0        0.7871   326        0.6792   ZeroHalf: 3     2037   15.90%
    *     0+    0                            0.6990        0.7871            12.61%
    *     0+    0                            0.7021        0.7871            12.12%
          0     2        0.7871   326        0.7021        0.7871     2037   12.12%
    Elapsed time = 0.62 sec. (488.67 ticks, tree = 0.00 MB, solutions = 7)
    *    11+    7                            0.7106        0.7636             7.46%
         60    41        0.7320   283        0.7106        0.7545     7871    6.18%
    *   135    83      integral     0        0.7158        0.7545    11716    5.42%
    *   137    83      integral     0        0.7160        0.7545    11763    5.38%
    *   172    80      integral     0        0.7232        0.7545    13516    4.33%
    
    Implied bound cuts applied:  55
    Zero-half cuts applied:  10
    Lift and project cuts applied:  1
    Gomory fractional cuts applied:  4
    
    Root node processing (before b&c):
      Real time             =    0.61 sec. (487.43 ticks)
    Parallel b&c, 8 threads:
      Real time             =    0.34 sec. (418.36 ticks)
      Sync time (average)   =    0.12 sec.
      Wait time (average)   =    0.13 sec.
                              ------------
    Total (root+branch&cut) =    0.95 sec. (905.80 ticks)
    
    Solution pool: 11 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  7.2320364849e-01
    Solution time =    0.95 sec.  Iterations = 31523  Nodes = 794
    Deterministic time = 905.80 ticks  (954.54 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpqt3_05gs.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.7232036484944747


We take the locations from the optimal monometallic problem to
initialize a **Canvas** object for the bimetallic case.

.. code:: ipython3

    Canv = Canvas()
    for i in range(len(D)):
        if(D.Contents[i] is not None):
            Canv.addLocation(D.Canvas.Points[i])
    Canv.setNeighborsFromFunc(Lat.getNeighbors)

Additionally, we create a few data structures for holding bimetallic
material information. First, we make a list of multiple **Atom** objects
that will be the building blocks of the model. Next, we specify a
dictionary with the bounds to impose on composition.

.. code:: ipython3

    Atoms = [Atom('Cu'),Atom('Ag')]
    CompBounds = {Atom('Cu'):(6,6),
                  Atom('Ag'):(14,14)}

Specifying an Optimization Model
--------------------------------

We start by creating a **MatOptModel** object that will hold the
information about the problem variables and constraints. At a minimum,
ever model requires a Canvas object to be defined. Additionally, the
list of building blocks and conformations that are present in the model
should be defined.

.. code:: ipython3

    m = MatOptModel(Canv,Atoms)

By default, several basic variables are pre-defined. See the first
example, **Monometallic_Nanocluster_Design.ipynb** for a description of
basic variables, expressions, and constraint rules.

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

    D = m.maximize(m.Ecoh,tilim=360,trelim=4096)


.. parsed-literal::

    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 12.6.1.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2014.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmp2mk8tqz0.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 360
    CPLEX> New value for upper limit on size of tree in megabytes: 4096
    CPLEX> Problem '/tmp/tmpbyvd2xuk.pyomo.lp' read.
    Read time = 0.00 sec. (0.07 ticks)
    CPLEX> Problem name         : /tmp/tmpbyvd2xuk.pyomo.lp
    Objective sense      : Maximize
    Variables            :     364  [Nneg: 1,  Fix: 2,  Free: 65,  Binary: 296]
    Objective nonzeros   :       1
    Linear constraints   :     876  [Less: 788,  Equal: 88]
      Nonzeros           :    2300
      RHS nonzeros       :     297
    
    Variables            : Min LB: 0.000000         Max UB: 14.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.01443376       Max   : 3.197034       
      RHS nonzeros       : Min   : 1.000000         Max   : 1.000000       
    CPLEX> Tried aggregator 2 times.
    MIP Presolve eliminated 87 rows and 68 columns.
    Aggregator did 20 substitutions.
    Reduced MIP has 769 rows, 276 columns, and 1812 nonzeros.
    Reduced MIP has 276 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (2.63 ticks)
    Found incumbent of value 1.960007 after 0.00 sec. (4.85 ticks)
    Probing time = 0.00 sec. (2.14 ticks)
    Tried aggregator 1 time.
    Reduced MIP has 769 rows, 276 columns, and 1812 nonzeros.
    Reduced MIP has 276 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (1.45 ticks)
    Probing time = 0.00 sec. (2.12 ticks)
    Clique table members: 1252.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 8 threads.
    Root relaxation solution time = 0.00 sec. (6.52 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
    *     0+    0                            1.9600        8.5512           336.28%
    *     0+    0                            2.0023        8.5512           327.07%
          0     0        3.3117   276        2.0023        3.3117      385   65.39%
    *     0+    0                            2.0587        3.3117            60.86%
          0     0        2.8575   210        2.0587     Cuts: 237      525   38.80%
          0     0        2.6095   176        2.0587      Cuts: 98      613   26.76%
          0     0        2.5198   150        2.0587      Cuts: 47      657   22.40%
          0     0        2.4525   208        2.0587      Cuts: 38      692   19.13%
          0     0        2.4119   150        2.0587      Cuts: 66      721   17.16%
          0     0        2.3480   174        2.0587      Cuts: 51      757   14.05%
          0     0        2.3266   136        2.0587      Cuts: 68      781   13.01%
          0     0        2.2652   169        2.0587      Cuts: 72      821   10.03%
          0     0        2.2146   148        2.0587      Cuts: 34      850    7.57%
    *     0+    0                            2.0815        2.2146             6.40%
          0     0        2.1706   117        2.0815      Cuts: 62      877    4.28%
          0     0        2.1405   139        2.0815      Cuts: 24      905    2.83%
          0     0        2.1276   171        2.0815      Cuts: 26      916    2.22%
          0     0        2.1106   185        2.0815      Cuts: 53      933    1.40%
          0     0        2.1079   160        2.0815      Cuts: 20      942    1.27%
    *     0+    0                            2.0871        2.1079             1.00%
          0     0        2.1032   160        2.0871      Cuts: 34      956    0.77%
          0     0        2.1016   182        2.0871  ZeroHalf: 17      976    0.70%
    *     0+    0                            2.0902        2.1016             0.55%
          0     0        2.0984   146        2.0902  ZeroHalf: 18      992    0.39%
    *     0+    0                            2.0958        2.0984             0.12%
          0     0        2.0967   160        2.0958  ZeroHalf: 23     1000    0.04%
          0     0        cutoff              2.0958                   1004    0.00%
    Elapsed time = 0.20 sec. (148.03 ticks, tree = 0.00 MB, solutions = 7)
    
    Clique cuts applied:  55
    Implied bound cuts applied:  62
    Zero-half cuts applied:  164
    
    Root node processing (before b&c):
      Real time             =    0.20 sec. (148.07 ticks)
    Parallel b&c, 8 threads:
      Real time             =    0.00 sec. (0.00 ticks)
      Sync time (average)   =    0.00 sec.
      Wait time (average)   =    0.00 sec.
                              ------------
    Total (root+branch&cut) =    0.20 sec. (148.07 ticks)
    
    Solution pool: 7 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  2.0958182296e+00
    Solution time =    0.20 sec.  Iterations = 1004  Nodes = 0
    Deterministic time = 148.07 ticks  (744.83 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpdsvgbv72.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 2.0958182295761083


Processing Solutions
--------------------

If a design was identified (optimal or otherwise), then a **Design**
object is returned from the optimization method. The optimal design can
be plotted via any of the supported parsers.

.. code:: ipython3

    if(D is not None):
        D.toPDB('result.pdb')
