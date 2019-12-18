Surface Design
==============

This notebook serves as an example application of the MatOpt framework.
We consider an example optimization problem of designing a monometallic
nanostructured catalyst surface.

Importing Packages
------------------

We start by importing several standard Python modules for convienience.

.. code:: ipython3

    import os 
    from math import sqrt
    import numpy as np
    from copy import deepcopy

Next, we import MatOpt.

.. code:: ipython3

    from matopt import *


::


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-2-be5adc2fb850> in <module>
    ----> 1 from matopt import *
    

    ModuleNotFoundError: No module named 'matopt'


Representing Materials
----------------------

To begin, we define a **Lattice** object. In this example,
**FCCLattice** is the appropriate a child class of Lattice. This object
will serve to define neighbor connections and helps us generically
create other objects. We construct our lattice from a class method
constructor for FCC lattices aligned with the {111} plane.

.. code:: ipython3

    IAD = 2.828427 # Angstrom
    Lat = FCCLattice.alignedWith111(IAD)

Next, we define a **Shape** object that we will use to specify a design
space. Additionally, in this example our design space is periodic, so we
will define a **Tiling** object to hold information about the
periodicity. In this example, **Parallelepiped** and **PlanarTiling**
are the appropriate child classes for these objects, respectively.

Note that we shift the shape of our design space slightly, in order to
avoid confusion about which lattice sites that lie perfectly on the
shape facet should be included.

.. code:: ipython3

    nUnitCellsOnEdge = 4
    nLayers = 6
    a = nUnitCellsOnEdge*IAD
    b = a
    c = nLayers*Lat.FCC111LayerSpacing
    alpha = np.pi/2
    beta = np.pi/2
    gamma = np.pi/3
    S = Parallelepiped.fromEdgesAndAngles(a,b,c,alpha,beta,gamma)
    S.shift(np.array([-0.01*a,-0.01*b,-0.01*c]))
    T = PlanarTiling(S)

Given the parameters for a design space, we can construct a **Canvas**
object to hold information about points and nearest neighbors. In this
example, the object is efficiently constructed from a scan over lattice
sites. In general, the Canvas can be constructed and manipulated via
user-defined algorithms.

.. code:: ipython3

    Canv = Canvas.fromLatticeAndTilingScan(Lat,T)

The Canvas object hold information about the design space and the
lattice sites, but it does not specify any material building block
information. To represent material configurations, use a **Design**
object.

Initially, the Design is empty. There are several ways to place **Atom**
(i.e., building block) objects in a Design. In this example, we are
considering a Pt surface for the oxygen reduction reaction. We can
initialize a design representing the FCC {111} surface by using a
standard constructor for the Design object.

To debug our work so far, we can create material structure files to load
and plot with standard visualization tools such as AtomEye. Here, we
create PDB (protein data bank format, www.rcsb.org) and CFG (AtomEye
configuration, li.mit.edu/A/Graphics/A/) files for the design. These
files can be plotted with visualization packages such as AtomEye or
OVITO.

.. code:: ipython3

    D = Design(Canv,Atom('Pt'))
    D.toPDB('undefected.pdb')
    D.toCFG('undefected.cfg',GS=1.0,BBox=S)

Building a Model
----------------

In this example, we will build a model that maximizes the number of
sites that are reactive for the oxygen reduction reaction (ORR). More
generally, our model will indicate sites that are within a certain
tolerance of a target generalized coordination number (GCN). These
target sites can also be constrained to lie within minimum and maximum
coordination number to be considered surface sites.

Additionally, we model the surface energy of nanostructured designs.
This surface energy can be constrained to be below a threshold and can
be included in the objective function. We can parametrically optimize
the multi-objective optimization problem by defining a weighting,
*CatWeight*, that controls how much weight is given to the catalytic
activity term in the objective function. A weighting of 1 corresponds to
the optimally active material and a weighting of 0 corresponds to the
lowest surface energy design.

.. code:: ipython3

    Atoms = [Atom('Pt')]
    TargetGCN = 8.0
    CNsurfMin = 3
    CNsurfMax = 9
    TileSizeSquared = nUnitCellsOnEdge**2
    UndefectedSurfE = 0.129758
    maxSurfE = 999
    CatWeight = 1.0

To begin, we start by creating a **MatOptModel** object to hold
information about the model.

.. code:: ipython3

    m = MatOptModel(Canv,Atoms)

By default, several basic variables are pre-defined. See the first
example, **Monometallic_Nanocluster_Design.ipynb** for a description of
basic variables, expressions, and constraint rules.

First, we introduce two rules to fix special sites in the design. We fix
the bottom two layers of atoms to exist, creating underlying bulk layers
above which we will introduce nanostruced defets. We also fix an
arbitrary atom in the top layer, breaking symetry of the design space
and resulting in easier to solve opitmization problems without actually
restricting the designs that can be possibly represented.

.. code:: ipython3

    CanvTwoBotLayers = [i for i in range(len(Canv)) 
                        if Canv.Points[i][2] < 1.5*Lat.FCC111LayerSpacing]
    CanvMinusTwoBotLayers = [i for i in range(len(Canv)) 
                             if i not in CanvTwoBotLayers]
    OneSiteInTopLayer = [min([i for i in range(len(Canv)) 
                              if Canv.Points[i][2] > (nLayers-1.5)*Lat.FCC111LayerSpacing])]
    m.Yi.rules.append(FixedTo(1,sites=OneSiteInTopLayer))
    m.Yi.rules.append(FixedTo(1,sites=CanvTwoBotLayers))

Next, we introduce constraints thtat require atoms to be placed on top
of each other, avoiding hollow pockets below the surface.

.. code:: ipython3

    NeighborsBelow = [[j for j in Canv.NeighborhoodIndexes[i] 
                       if(j is not None and
                          Canv.Points[j][2]<Canv.Points[i][2]-DBL_TOL)] 
                      for i in range(len(Canv))]
    m.Yi.rules.append(ImpliesNeighbors(concs=(m.Yi,GreaterThan(1)),
                                       sites=CanvMinusTwoBotLayers,
                                       neighborhoods=NeighborsBelow))

Next, we introduce several rules for the geometric and reactive
descriptors of sites in the design. We define the generalized
coordination number according to a linear equality constraint. Then, we
define ideal sites as having a conjunction of requirements on the
generalized coordination number, and regular coordination number.
Finally, we define activity as the count of sites with target
coordination number.

.. code:: ipython3

    m.addSitesDescriptor('GCNi',bounds=(0,12),integer=False,
                         rules=EqualTo(SumNeighborSites(desc=m.Ci,
                                                        coefs=1/12)),
                         sites=CanvMinusTwoBotLayers)
    m.addSitesDescriptor('IdealSitei',binary=True,
                         rules=[Implies(concs=(m.Ci,GreaterThan(3))),
                                Implies(concs=(m.Ci,LessThan(9))),
                                Implies(concs=(m.GCNi,EqualTo(TargetGCN)))],
                         sites=CanvMinusTwoBotLayers)
    m.addGlobalDescriptor('Activity',bounds=(0,1),
                          rules=EqualTo(SumSites(m.IdealSitei,coefs=1/TileSizeSquared)))

Next, we define a simple model for the surface energy of nanostructured
slabs as a piecwise linear function of coordination number.

.. code:: ipython3

    EiVals = [0, -0.04293*3+0.41492, -0.04293*10+0.41492, 0.05179*11-0.62148, 0]
    EiBPs = [0, 3, 10, 11, 12]
    m.addSitesDescriptor('Ei',rules=PiecewiseLinear(values=EiVals,
                                                    breakpoints=EiBPs,
                                                    input_desc=m.Ci),
                         sites=CanvMinusTwoBotLayers)
    m.addGlobalDescriptor('Esurf',bounds=(None,maxSurfE),
                          rules=EqualTo(SumSites(m.Ei,coefs=1/TileSizeSquared,offset=0.101208)))
    m.addGlobalDescriptor('Stability',
                          rules=EqualTo(LinearExpr(m.Esurf,1/UndefectedSurfE)))

Finally, we introduce a descriptor for the weighted combination of
activity and stability.

.. code:: ipython3

    m.addGlobalDescriptor('ActAndStab',
                          rules=EqualTo(LinearExpr(descs=[m.Stability,m.Activity],
                                                          coefs=[-(1-CatWeight),CatWeight])))

Solving the Model
-----------------

Given a fully formed Pyomo model, we have several capabilities to
optimize and visualize the solution. In this example, we simply call the
maximize method to optimize the balance of activity and stability

.. code:: ipython3

    D = m.maximize(m.ActAndStab,tilim=360)


.. parsed-literal::

    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 12.6.1.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2014.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmppgmxz9vy.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 360
    CPLEX> Problem '/tmp/tmpmqh62f6a.pyomo.lp' read.
    Read time = 0.01 sec. (0.32 ticks)
    CPLEX> Problem name         : /tmp/tmpmqh62f6a.pyomo.lp
    Objective sense      : Maximize
    Variables            :    1428  [Nneg: 321,  Box: 65,  Free: 66,  Binary: 895,
                                     General Integer: 80,  Other: 1]
    Objective nonzeros   :       1
    Linear constraints   :    3573  [Less: 3232,  Equal: 341]
      Nonzeros           :    8656
      RHS nonzeros       :    1029
    SOS                  :      64  [SOS2: 64, 320 members, all continuous]
    
    Variables            : Min LB: 0.000000         Max UB: 999.0000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.01438000       Max   : 12.00000       
      RHS nonzeros       : Min   : 0.1012080        Max   : 12.00000       
    CPLEX> MIP Presolve eliminated 1 redundant SOS constraints.
    Tried aggregator 2 times.
    MIP Presolve eliminated 2237 rows and 539 columns.
    MIP Presolve modified 488 coefficients.
    Aggregator did 256 substitutions.
    Reduced MIP has 1080 rows, 633 columns, and 4051 nonzeros.
    Reduced MIP has 301 binaries, 66 generals, 63 SOSs, and 0 indicators.
    Presolve time = 0.01 sec. (10.00 ticks)
    Probing fixed 0 vars, tightened 70 bounds.
    Probing time = 0.03 sec. (30.84 ticks)
    Tried aggregator 1 time.
    Reduced MIP has 1080 rows, 633 columns, and 4051 nonzeros.
    Reduced MIP has 301 binaries, 66 generals, 63 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (2.83 ticks)
    Probing fixed 0 vars, tightened 6 bounds.
    Probing time = 0.02 sec. (23.29 ticks)
    Clique table members: 2849.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 8 threads.
    Root relaxation solution time = 0.01 sec. (13.74 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
          0     0        1.0000   182                      1.0000      196         
    *     0+    0                            0.0000        1.0000              --- 
          0     0        1.0000   148        0.0000      Cuts: 13      249     --- 
          0     0        1.0000   269        0.0000     Cuts: 167      427     --- 
    *     0+    0                            0.0625        1.0000              --- 
          0     0        1.0000   187        0.0625      Cuts: 15      460     --- 
          0     0        1.0000   310        0.0625      Cuts: 96      620     --- 
          0     2        1.0000   274        0.0625        1.0000      620     --- 
    Elapsed time = 0.47 sec. (420.18 ticks, tree = 0.00 MB, solutions = 2)
    *   281+  116                            0.1875        1.0000           433.33%
        296    60        0.8914   209        0.1875        1.0000    10764  433.33%
    *   478+  123                            0.2500        1.0000           300.00%
        770   161        0.8694   239        0.2500        1.0000    30637  300.00%
       1219   303        0.9495   289        0.2500        1.0000    57291  300.00%
    *  1580   432      integral     0        0.3750        0.9460    73750  152.27%
       1833   217        0.4420   128        0.3750        0.8149    88172  117.32%
    
    Clique cuts applied:  6
    Cover cuts applied:  17
    Implied bound cuts applied:  249
    Mixed integer rounding cuts applied:  1
    Zero-half cuts applied:  3
    
    Root node processing (before b&c):
      Real time             =    0.47 sec. (419.13 ticks)
    Parallel b&c, 8 threads:
      Real time             =    1.16 sec. (1061.27 ticks)
      Sync time (average)   =    0.27 sec.
      Wait time (average)   =    0.28 sec.
                              ------------
    Total (root+branch&cut) =    1.63 sec. (1480.40 ticks)
    
    Solution pool: 6 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  3.7500000000e-01
    Solution time =    1.63 sec.  Iterations = 98707  Nodes = 2303
    Deterministic time = 1480.40 ticks  (908.51 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpqdq58zzx.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.375


Processing Solutions
--------------------

Once the model is solved, we can plot the resulting design. However, it
is often useful to label atoms according to some auxilliary information.
In this case, we would like to label atoms that consitute ideal reactive
sites. We loop over the sites and set the atom to S to highlight the
sites that are reactive. Then, we can write the Design object to PDB or
CFG files for plotting.

Additionally, we can manipulate the resulting design to better see the
periodic pattern. Here, we replicate the design four times to see the
periodic pattern.

.. code:: ipython3

    if(D is not None):
        for i in m.IdealSitei.keys():
            if m.IdealSitei.values[i] > 0.5:
                D.setContent(i,Atom('S'))
        D.toPDB('result.pdb')
        PeriodicD = T.replicateDesign(D,4)
        PeriodicS = deepcopy(S)
        PeriodicS.scale(np.array([4,4,1]))
        PeriodicD.toCFG('periodic_result.cfg',BBox=PeriodicS)
