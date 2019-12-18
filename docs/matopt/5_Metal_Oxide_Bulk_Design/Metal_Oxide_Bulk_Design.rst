Perovskite Design
=================

This notebook serves as an example application of the MatOpt framework
on bulk nanostructured materials. We consider the problem of how to
optimally place dopant in a perovskite lattice.

For more information, see: Hanselman, Christopher L., et al. “A
framework for optimizing oxygen vacancy formation in doped perovskites.”
*Computers & Chemical Engineering* 126 (2019): 168-177. DOI:
`10.1016/j.compchemeng.2019.03.033 <https://doi.org/10.1016/j.compchemeng.2019.03.033>`__

Importing Packages
------------------

We start by importing MatOpt.

.. code:: ipython3

    from matopt import *


::


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-1-be5adc2fb850> in <module>
    ----> 1 from matopt import *
    

    ModuleNotFoundError: No module named 'matopt'


Representing Materials
----------------------

First, we construct a **Lattice** object to hold information about the
sites in our material. For this application, we use the
**PerovskiteLattice** class with lattice constants *A*, *B*, and *C*.

.. code:: ipython3

    A = 4.0
    B = 4.0
    C = 4.0
    Lat = PerovskiteLattice(A,B,C)

Next, we construct **Shape** and **Tiling** objects to help define the
material locations of interest. In this case, we use **RectPrism** and
**CubicTiling**, respectively.

Note that we shift the shape slightly to avoid ambiguity about which
sites on the border of the cell should be included.

.. code:: ipython3

    nUnitCellsOnEdge = 2
    S = RectPrism(nUnitCellsOnEdge*A,
                  nUnitCellsOnEdge*B,
                  nUnitCellsOnEdge*C)
    S.shift(np.array([-0.01,-0.01,-0.01]))
    T = CubicTiling(S)

Next, we construct the **Canvas** object which will hold all the
information about the sites and information about neighbors. We also
define a list of **Atom** objects that serve as the building blocks of
our material.

.. code:: ipython3

    Canv = Canvas.fromLatticeAndTilingScan(Lat,T)
    Atoms = [Atom('Ba'),Atom('Fe'),Atom('In'),Atom('O')]

Finally, we load a list of conformations from file that represent a set
of dopant configurations that we would like to indicate in the design.

.. code:: ipython3

    iDesiredConfs = [394,395,396,397,398,399,400,401,68,69,
                     70,71,162,163,164,165,166,167,168,169]
    ConfDesigns = loadFromPDBs([str(i)+'.pdb' for i in iDesiredConfs],folder='./Confs')
    Confs = [Conf.Contents for Conf in ConfDesigns]

Building the Model
------------------

To begin specifying the model, we first define several pieces of
information that will help specify the design problem.

.. code:: ipython3

    Sites = [i for i in range(len(Canv))]
    ASites = [i for i in Sites if Lat.isASite(Canv.Points[i])]
    BSites = [i for i in Sites if Lat.isBSite(Canv.Points[i])]
    OSites = [i for i in Sites if Lat.isOSite(Canv.Points[i])]
    pctLocalLB,pctLocalUB = 0,1
    pctGlobalLB,pctGlobalUB = 0.0,0.3
    LocalBounds = {(i,Atom('In')):(round(pctLocalLB*len(Canv.NeighborhoodIndexes[i])),
                                   round(pctLocalUB*len(Canv.NeighborhoodIndexes[i]))) for i in OSites}
    GlobalLB = round(pctGlobalLB*len(BSites))
    GlobalUB = round(pctGlobalUB*len(BSites))

Next, we initialize a **MatOptModel** object that will hold all the
information about material descriptors and desired functionalities.

.. code:: ipython3

    m = MatOptModel(Canv,Atoms,Confs)

By default, several basic variables are pre-defined. See the first
example, **Monometallic_Nanocluster_Design.ipynb** for a description of
basic variables, expressions, and constraint rules.

For this system, we introduce several rules about the allowed placement
of atoms in the design. First, we require that all A-sites in the
material are occupied by Ba. Next, we require that all O-sites are
occupied by O. Thirdly, we forbid Ba and O from being placed in B-sites.
And finally, we require that some atom be placed in each B-site. These
four rules effectively limit the scope of the optimization to focus on
the labeling of B-sites as either Fe or In.

.. code:: ipython3

    m.Yik.rules.append(FixedTo(1,sites=ASites,site_types=[Atom('Ba')]))
    m.Yik.rules.append(FixedTo(1,sites=OSites,site_types=[Atom('O')]))
    m.Yik.rules.append(FixedTo(0,sites=BSites,site_types=[Atom('Ba'),Atom('O')]))
    m.Yi.rules.append(FixedTo(1,sites=BSites))

To specify additional constraints to the model, we create several
descriptors for the activity, local dopant concentration, and the global
dopant concentration.

Notice that in each case, we specify a subset of locations or atoms of
interest. This is because, for example, our material activity depends on
oxygen sites only and it would be nonsensical to try to interpret one of
the conformations on a different type of site. Similarly, the dopant
budgets are written only over In atoms and not on Ba, Fe, or O.

.. code:: ipython3

    m.addGlobalDescriptor('Activity',
                          rules=EqualTo(SumSitesAndConfs(m.Zic,coefs=1/len(OSites),sites_to_sum=OSites)))
    m.addGlobalTypesDescriptor('GlobalIndiumConc',bounds=(GlobalLB,GlobalUB),
                               rules=EqualTo(SumSites(m.Yik,
                                                      site_types=[Atom('In')],
                                                      sites_to_sum=BSites)))
    m.addSitesTypesDescriptor('LocalIndiumConc',bounds=LocalBounds,
                              rules=EqualTo(SumNeighborSites(m.Yik,
                                                             sites=OSites,
                                                             site_types=[Atom('In')])))

Solving the Model
-----------------

Given a fully formed model, we can optimize by maximizing or minimizing
one of the global descriptors.

.. code:: ipython3

    D = m.maximize(m.Activity,tilim=360)


.. parsed-literal::

    
    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 12.6.1.0
      with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2014.  All Rights Reserved.
    
    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.
    
    CPLEX> Logfile 'cplex.log' closed.
    Logfile '/tmp/tmpp1ouabbf.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 360
    CPLEX> Problem '/tmp/tmpv8p2na6w.pyomo.lp' read.
    Read time = 0.01 sec. (0.32 ticks)
    CPLEX> Problem name         : /tmp/tmpv8p2na6w.pyomo.lp
    Objective sense      : Maximize
    Variables            :     843  [Nneg: 1,  Box: 25,  Free: 1,  Binary: 816]
    Objective nonzeros   :       1
    Linear constraints   :    1459  [Less: 1424,  Equal: 35]
      Nonzeros           :   17523
      RHS nonzeros       :    1113
    
    Variables            : Min LB: 0.000000         Max UB: 10.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.04166667       Max   : 2.000000       
      RHS nonzeros       : Min   : 1.000000         Max   : 9.000000       
    CPLEX> Found incumbent of value 0.000000 after 0.00 sec. (0.18 ticks)
    Tried aggregator 2 times.
    MIP Presolve eliminated 1042 rows and 515 columns.
    MIP Presolve modified 3072 coefficients.
    Aggregator did 8 substitutions.
    Reduced MIP has 409 rows, 320 columns, and 3152 nonzeros.
    Reduced MIP has 320 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.01 sec. (14.68 ticks)
    Probing fixed 144 vars, tightened 0 bounds.
    Probing time = 0.01 sec. (7.84 ticks)
    Tried aggregator 1 time.
    MIP Presolve eliminated 3 rows and 147 columns.
    Reduced MIP has 406 rows, 173 columns, and 2123 nonzeros.
    Reduced MIP has 173 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (2.75 ticks)
    Probing time = 0.00 sec. (4.49 ticks)
    Tried aggregator 1 time.
    MIP Presolve eliminated 2 rows and 2 columns.
    Reduced MIP has 404 rows, 171 columns, and 2109 nonzeros.
    Reduced MIP has 171 binaries, 0 generals, 0 SOSs, and 0 indicators.
    Presolve time = 0.00 sec. (2.61 ticks)
    Probing time = 0.00 sec. (3.09 ticks)
    Clique table members: 955.
    MIP emphasis: balance optimality and feasibility.
    MIP search method: dynamic search.
    Parallel mode: deterministic, using up to 8 threads.
    Root relaxation solution time = 0.00 sec. (3.99 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
    *     0+    0                            0.0000        7.0000              --- 
    *     0+    0                            0.3333        7.0000              --- 
          0     0        1.0000   101        0.3333        1.0000      326  200.00%
    *     0+    0                            0.5000        1.0000           100.00%
          0     0        0.9714   118        0.5000     Cuts: 168      452   94.27%
          0     0        0.9167   106        0.5000     Cuts: 159      522   83.33%
          0     0        0.9033   106        0.5000      Cuts: 38      573   80.65%
          0     0        0.8907   124        0.5000      Cuts: 58      672   78.15%
          0     0        0.8686   130        0.5000  ZeroHalf: 65      773   73.72%
          0     0        cutoff              0.5000        0.5000      773    0.00%
    Elapsed time = 0.10 sec. (96.02 ticks, tree = 0.00 MB, solutions = 3)
    
    Clique cuts applied:  13
    Zero-half cuts applied:  31
    Lift and project cuts applied:  1
    Gomory fractional cuts applied:  4
    
    Root node processing (before b&c):
      Real time             =    0.10 sec. (96.07 ticks)
    Parallel b&c, 8 threads:
      Real time             =    0.00 sec. (0.00 ticks)
      Sync time (average)   =    0.00 sec.
      Wait time (average)   =    0.00 sec.
                              ------------
    Total (root+branch&cut) =    0.10 sec. (96.07 ticks)
    
    Solution pool: 4 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  5.0000000000e-01
    Solution time =    0.10 sec.  Iterations = 773  Nodes = 0
    Deterministic time = 96.07 ticks  (981.48 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmpckxysoyl.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.4999999999999998


Processing Solutions
--------------------

If the optimizer was successful in finding an optimal (or just feasible)
solution, we can plot the resulting design to any of several standard
file formats. However, it is often useful to modify the design to
highlight key features. Here, we label all O-sites that constitute one
of the desired conformations by replacing the atom with an S.

.. code:: ipython3

    if(D is not None):
        for i,c in m.Zic.keys():
            if(m.Zic.values[i,c] > 0.5):
                D.setContent(i,Atom('S'))
        D.toCFG('result.cfg',BBox=S)
