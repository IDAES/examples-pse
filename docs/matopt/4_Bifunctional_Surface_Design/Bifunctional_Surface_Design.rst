Bifunctional Catalyst Design
============================

This notebook serves as an example application of the MatOpt framework.
We consider an example optimization problem of designing a
nanostructured bifunctional catalyst. This example is a simplified
representation of the system presented in [1].

[1] Nunez, M., & Vlachos, D. G. (2019). Ind. Eng. Chem. Res., 58,
6146-6154.

Importing Packages
------------------

We start by importing several standard Python modules for convienience.

.. code:: ipython3

    import os 
    from copy import deepcopy
    import numpy as np

Finally, we import the MatOpt package in its entirety.

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
**FCCLattice** is a child class of Lattice. This object will serve to
define neighbor connections and helps us generically create other
objects.

.. code:: ipython3

    IAD = 2.828 # Angstrom
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

    nUnitCellsOnEdge = 8
    nLayers = 4
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
initialize the Design to hold all Pt atoms.

To debug our work so far, we can create material structure files to load
and plot with standard visualization tools such as AtomEye. Here, we
create PDB (protein data bank format, www.rcsb.org) and CFG (AtomEye
configuration, li.mit.edu/A/Graphics/A/) files for the undoped design.

.. code:: ipython3

    D = Design(Canv,Atom('Pt'))
    D.toPDB('canvas.pdb')
    D.toCFG('canvas.cfg',BBox=S)

Representing Conformations
--------------------------

In this material system, we would like to model the presence of facet
and edge sites on a patchy bimetallic catalyst surface. To do this
generically, we will create a list of conformations. This list will
later be used by MatOpt modeling methods to create common descriptor
formulations.

To begin, we create another Canvas object with one shell of neighbors
around a lattice location. Then, we create a list of Designs and set
their contents to match our intended conformations. To debug our work,
we also output conformations to file for plotting.

.. code:: ipython3

    MotifCanvas = Canvas()
    MotifCanvas.addLocation(np.array([0,0,0],dtype=float),NNeighbors=12)
    MotifCanvas.addShell(Lat.getNeighbors)
    Confs = [[None]*len(MotifCanvas.NeighborhoodIndexes[0]) for _ in range(7)]
    iToSetNi = [[3,4,5,6,7,8],
                [3,4,5,6],
                [4,5,6,7],
                [5,6,7,8],
                [6,7,8,3],
                [7,8,3,4],
                [8,3,4,5]]
    iToSetPt = [[9,10,11],
                [9,10,11],
                [9,10,11],
                [9,10,11],
                [9,10,11],
                [9,10,11],
                [9,10,11]]
    for iConf,Conf in enumerate(Confs):
        for i in iToSetNi[iConf]:
            Conf[i] = Atom('Ni')
        for i in iToSetPt[iConf]:
            Conf[i] = Atom('Pt')


Building the Model
------------------

To begin, we define several sets and constants that will be used in
creating the model.

.. code:: ipython3

    TypeAConfs = [0]
    TypeBConfs = [1,2,3,4,5,6]
    LocsToFixPt = [i for i in range(len(Canv)) if Canv.Points[i][2] < Lat.FCC111LayerSpacing*2.5]
    LocsToExcludePt = [i for i in range(len(Canv)) if i not in LocsToFixPt]
    CanvTwoBotLayers = [i for i in range(len(Canv)) if Canv.Points[i][2] < Lat.FCC111LayerSpacing*1.5]
    CanvMinusTwoBotLayers = [i for i in range(len(Canv)) if i not in CanvTwoBotLayers]
    OneLocToFix = [min(LocsToExcludePt)]
    TileSizeSquared = nUnitCellsOnEdge**2
    CatNorm = TileSizeSquared*6.0
    UndefectedSurfE = 0.129758
    maxSurfE = 999
    CatWeight = 1.0
    Atoms = [Atom('Ni'),Atom('Pt')]

Next, we create a **MatOptModel** object.

.. code:: ipython3

    m = MatOptModel(Canv,Atoms,Confs)

By default, several basic variables are pre-defined. See the first
example, **Monometallic_Nanocluster_Design.ipynb** for a description of
basic variables, expressions, and constraint rules.

First, we fix the composition of atoms in the appropriate layers.
Effectively, we are designing the defects in a single layer of Ni on top
of an undefected Pt surface.

.. code:: ipython3

    m.Yik.rules.append(FixedTo(1,sites=LocsToFixPt,site_types=[Atom('Pt')]))
    m.Yik.rules.append(FixedTo(0,sites=LocsToExcludePt,site_types=[Atom('Pt')]))

Next, we define indicators for the presence of groups of conformations
(corresponding to facet and edge sites) in the design. We arbitrarily
fix one site to be a facet-type site, breaking symmetry and improving
the tractability of the resulting optimization models.

.. code:: ipython3

    m.Zic.rules.append(FixedTo(1,sites=OneLocToFix,confs=TypeAConfs))
    m.Zic.rules.append(Implies(concs=(m.Yik,EqualTo(1,site_types=[Atom('Ni')]))))
    SumAConfsExpr = SumConfs(m.Zic,confs_to_sum=TypeAConfs)
    SumBConfsExpr = SumConfs(m.Zic,confs_to_sum=TypeBConfs)
    m.addBondsDescriptor('SiteCombinations',binary=True,
                         rules=ImpliesSiteCombination(Canv,
                                                      (SumAConfsExpr,GreaterThan(1)),
                                                      (SumBConfsExpr,GreaterThan(1))))

Next, we define activity as a normalized sum of contributions from site
combinations. Additionally, we introduce a model for the surface energy
of sites as a piecewise linear function of coordination number.

.. code:: ipython3

    m.addGlobalDescriptor('Activity',
                          rules=EqualTo(SumBonds(m.SiteCombinations,coefs=1/CatNorm)))
    
    EiVals = [0, -0.04293*3+0.41492, -0.04293*10+0.41492, 0.05179*11-0.62148, 0]
    EiBPs = [0, 3, 10, 11, 12]
    m.addSitesDescriptor('Ei',
                         rules=PiecewiseLinear(values=EiVals,
                                               breakpoints=EiBPs,
                                              input_desc=m.Ci),
                         sites=CanvMinusTwoBotLayers)
    m.addGlobalDescriptor('Esurf',
                          rules=EqualTo(SumSites(m.Ei,coefs=1/TileSizeSquared,offset=0.101208)))
    m.addGlobalDescriptor('Stability',
                          rules=EqualTo(LinearExpr(m.Esurf,1/UndefectedSurfE)))

Finally, we introduce a single descriptor for the weighted combination
of acitivity and stability. By changing the parameter weighting the
catalytic portion of the objective function, we can optimize for a range
of designs optimizing stability and activity.

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
    Logfile '/tmp/tmpe4o62som.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 360
    CPLEX> Problem '/tmp/tmplv6jodwi.pyomo.lp' read.
    Read time = 0.06 sec. (3.54 ticks)
    CPLEX> Problem name         : /tmp/tmplv6jodwi.pyomo.lp
    Objective sense      : Maximize
    Variables            :    7810  [Nneg: 641,  Free: 132,  Binary: 6653,
                                     General Integer: 384]
    Objective nonzeros   :       1
    Linear constraints   :   35589  [Less: 34304,  Equal: 1285]
      Nonzeros           :   90187
      RHS nonzeros       :   12680
    SOS                  :     128  [SOS2: 128, 640 members, all continuous]
    
    Variables            : Min LB: 0.000000         Max UB: 12.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.002604167      Max   : 12.00000       
      RHS nonzeros       : Min   : 0.1012080        Max   : 9.000000       
    CPLEX> MIP Presolve eliminated 10 redundant SOS constraints.
    Tried aggregator 3 times.
    MIP Presolve eliminated 31737 rows and 5492 columns.
    MIP Presolve modified 362 coefficients.
    Aggregator did 703 substitutions.
    Reduced MIP has 3149 rows, 1615 columns, and 12613 nonzeros.
    Reduced MIP has 1146 binaries, 57 generals, 118 SOSs, and 0 indicators.
    Presolve time = 0.03 sec. (40.22 ticks)
    Probing fixed 36 vars, tightened 99 bounds.
    Probing time = 0.03 sec. (14.93 ticks)
    Tried aggregator 1 time.
    MIP Presolve eliminated 72 rows and 36 columns.
    Reduced MIP has 3077 rows, 1579 columns, and 12409 nonzeros.
    Reduced MIP has 1110 binaries, 57 generals, 118 SOSs, and 0 indicators.
    Presolve time = 0.01 sec. (10.30 ticks)
    Probing fixed 0 vars, tightened 12 bounds.
    Probing time = 0.02 sec. (12.53 ticks)
    Clique table members: 6970.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 8 threads.
    Root relaxation solution time = 0.19 sec. (205.97 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
          0     0        0.3585  1118                      0.3585     2783         
    *     0+    0                            0.0000        0.3585              --- 
          0     0        0.3377  1101        0.0000     Cuts: 645     3852     --- 
          0     0        0.3101  1205        0.0000     Cuts: 785     5444     --- 
          0     0        0.2928  1254        0.0000     Cuts: 681     6792     --- 
          0     0        0.2844  1244        0.0000     Cuts: 399     7665     --- 
    *     0+    0                            0.0599        0.2844           374.86%
          0     0        0.2763  1229        0.0599     Cuts: 287     8590  361.36%
          0     0        0.2697  1228        0.0599     Cuts: 327     9467  350.21%
          0     0        0.2629  1209        0.0599     Cuts: 303    10404  338.98%
    *     0+    0                            0.0729        0.2629           260.59%
          0     0        0.2596  1216        0.0729     Cuts: 215    11274  256.07%
    *     0+    0                            0.0833        0.2596           211.56%
          0     0        0.2565  1211        0.0833     Cuts: 260    12013  207.83%
          0     0        0.2541  1215        0.0833     Cuts: 269    12740  204.94%
          0     0        0.2499  1208        0.0833     Cuts: 231    13818  199.94%
          0     0        0.2473  1200        0.0833     Cuts: 175    14508  196.79%
          0     0        0.2461  1210        0.0833     Cuts: 234    15153  195.26%
          0     0        0.2452  1206        0.0833     Cuts: 128    15706  194.20%
          0     0        0.2446  1204        0.0833     Cuts: 103    16244  193.49%
          0     0        0.2438  1209        0.0833     Cuts: 262    16892  192.60%
          0     0        0.2432  1217        0.0833     Cuts: 134    17461  191.85%
          0     0        0.2428  1213        0.0833     Cuts: 203    17924  191.31%
          0     0        0.2424  1218        0.0833     Cuts: 126    18473  190.88%
          0     0        0.2422  1223        0.0833      Cuts: 83    18844  190.61%
          0     0        0.2420  1212        0.0833     Cuts: 154    19227  190.39%
          0     0        0.2418  1222        0.0833     Cuts: 108    19629  190.16%
          0     0        0.2414  1215        0.0833     Cuts: 174    20328  189.64%
    *     0+    0                            0.0859        0.2414           180.86%
          0     0        0.2410  1217        0.0859     Cuts: 154    21041  180.45%
    *     0+    0                            0.0911        0.2410           164.42%
          0     0        0.2406  1214        0.0911     Cuts: 112    21706  163.94%
    *     0+    0                            0.0938        0.2406           156.61%
          0     0        0.2403  1219        0.0938     Cuts: 125    22279  156.32%
          0     0        0.2398  1214        0.0938     Cuts: 184    22915  155.81%
          0     0        0.2394  1215        0.0938     Cuts: 131    23512  155.40%
          0     0        0.2390  1220        0.0938     Cuts: 147    23976  154.96%
          0     0        0.2386  1223        0.0938     Cuts: 156    24511  154.48%
          0     0        0.2383  1217        0.0938      Cuts: 92    24950  154.17%
          0     0        0.2377  1214        0.0938     Cuts: 157    25548  153.58%
          0     0        0.2373  1215        0.0938     Cuts: 117    26062  153.11%
          0     0        0.2369  1213        0.0938     Cuts: 103    26470  152.73%
          0     0        0.2364  1217        0.0938     Cuts: 106    26964  152.19%
          0     0        0.2360  1212        0.0938     Cuts: 107    27523  151.74%
          0     0        0.2358  1207        0.0938      Cuts: 75    27893  151.52%
    *     0+    0                            0.1354        0.2358            74.13%
          0     0        0.2354  1207        0.1354      Cuts: 89    28464   73.84%
          0     0        0.2349  1208        0.1354     Cuts: 303    29117   73.45%
          0     0        0.2345  1215        0.1354     Cuts: 269    29791   73.16%
          0     0        0.2341  1213        0.1354      Cuts: 88    30279   72.90%
          0     0        0.2340  1215        0.1354     Cuts: 137    30662   72.80%
          0     2        0.2340  1204        0.1354        0.2320    30662   71.33%
    Elapsed time = 15.81 sec. (13534.85 ticks, tree = 0.00 MB, solutions = 8)
          1     3        0.2272  1105        0.1354        0.2320    32190   71.33%
          3     5        0.2326  1150        0.1354        0.2320    33587   71.33%
          6     8        0.2215  1084        0.1354        0.2320    36638   71.33%
          8    10        0.2179  1042        0.1354        0.2320    40007   71.33%
         10    12        0.2289  1149        0.1354        0.2320    44768   71.33%
         16    18        0.2178  1085        0.1354        0.2320    61328   71.33%
         22    24        0.2050   997        0.1354        0.2320    73062   71.33%
         37    39        0.1724   844        0.1354        0.2320    95678   71.33%
         62    60        0.2256  1116        0.1354        0.2320   116094   71.33%
        135   109        0.2075  1057        0.1354        0.2320   168727   71.33%
    Elapsed time = 19.92 sec. (17706.48 ticks, tree = 0.00 MB, solutions = 8)
        210   162        0.1615   837        0.1354        0.2320   225484   71.33%
    *   320+  197                            0.1667        0.2320            39.21%
        324   200        0.2021  1005        0.1667        0.2320   288547   39.21%
        383   135        0.1927  1004        0.1667        0.2273   344865   36.39%
        452   150        0.2211  1105        0.1667        0.2227   395287   33.61%
        518   182        cutoff              0.1667        0.2211   457859   32.67%
        548   200        0.1752   942        0.1667        0.2211   503022   32.67%
        600   222        0.1948  1002        0.1667        0.2211   554457   32.67%
        728   278        0.1867   977        0.1667        0.2211   592821   32.67%
        854   306        0.1753   957        0.1667        0.2034   658469   22.03%
        894   314        0.1885  1013        0.1667        0.2034   691508   22.03%
    Elapsed time = 32.31 sec. (27875.03 ticks, tree = 0.17 MB, solutions = 9)
       1016   338        0.1699   892        0.1667        0.1922   772028   15.34%
       1106   354        0.1733   909        0.1667        0.1922   813306   15.34%
       1200   346        0.1695   936        0.1667        0.1884   866432   13.06%
       1269   311        cutoff              0.1667        0.1884   923118   13.06%
       1381   244        0.1700   903        0.1667        0.1884  1011156   13.06%
       1494   169        cutoff              0.1667        0.1809  1067652    8.55%
    
    Clique cuts applied:  62
    Implied bound cuts applied:  1298
    Flow cuts applied:  1
    Mixed integer rounding cuts applied:  22
    Zero-half cuts applied:  280
    Lift and project cuts applied:  1
    
    Root node processing (before b&c):
      Real time             =   15.78 sec. (13503.33 ticks)
    Parallel b&c, 8 threads:
      Real time             =   24.60 sec. (20904.59 ticks)
      Sync time (average)   =    2.60 sec.
      Wait time (average)   =    2.62 sec.
                              ------------
    Total (root+branch&cut) =   40.37 sec. (34407.92 ticks)
    
    Solution pool: 9 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  1.6666666667e-01
    Solution time =   40.38 sec.  Iterations = 1096457  Nodes = 1666
    Deterministic time = 34407.94 ticks  (852.20 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpe4m74p2d.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.1666666666666666


Processing Solutions
--------------------

Once the model is solved, we can interpret the solutions as labelings of
a Design object. To accompolish this, we use the **setDesignFromModel**
function. Then, we can write the Design object to PDB or CFG files for
plotting.

.. code:: ipython3

    if(D is not None):
            D.toCFG('result.cfg',BBox=S)
            PeriodicD = T.replicateDesign(D,4)
            PeriodicS = deepcopy(S)
            PeriodicS.scale(np.array([4,4,1]))
            PeriodicD.toCFG('periodic_result.cfg',BBox=PeriodicS)
