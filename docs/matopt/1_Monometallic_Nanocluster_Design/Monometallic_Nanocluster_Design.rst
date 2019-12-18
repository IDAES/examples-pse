Monometallic Nanocluster Design
===============================

In this module, we introduce the **MatOpt** interface for representing
material properties and specifying optimization problems.

We have designed the interface with severl goals in mind:

1. To **simplify the representation of nanostructured materials,**
   streamlining the creation of materials optimization problems.
2. To provide a simple interface so that users **do not need to
   understand the details of building mathematical optmization models**
   or the syntax of the Pyomo package.
3. To **automate many of the common steps of materials optimization,**
   speeding up the development of new models.

As an example system, we will consider the minimization of cohesive
energy in nanoclusters, recently demonstrated in:

Isenberg, N. M., et al., “Identification of Optimally Stable Nanocluster
Geometries via Mathematical Optimization and Density Functional Theory,”
*Molecular Systems Design & Engineering*, 2019. DOI:
`10.1039/C9ME00108E <https://doi.org/10.1039/C9ME00108E>`__.

We seek to identify the geometry that minimizes the cohesive energy of a
nanocluster on the face-centered cubic (FCC) lattice. As a model for
cohesive energy, we use model based on the square-root of coordination
number, refered to as the Tomanek model
`[1] <https://doi.org/10.1103/PhysRevB.28.665>`__. In the equation
below, we define the normalized cohesive energy, as the normalized
contribution of the square root of the coordination number.

.. math:: \hat{E}^{\text{surf}} = \frac{1}{N \sqrt{12}} \displaystyle \sum_i \sqrt{CN_i} 

In the following sections, we demonstrate the basic approach for
importing the MatOpt package, specifying the design space, formulating
an optimization model, solving the optimization problem, and then
outputting results.

Importing Packages
------------------

We start by importing several standard Python modules for convienience.

.. code:: ipython3

    import os
    from math import sqrt

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

To begin formulating a material optimization problem, we need several
pieces of information about the design space. Our goal is to generate a
data structure for representing the choices in the design space, namely
the choice of where to place atoms on FCC lattice sites.

First, we define an **FCCLattice** object that holds the information
about what sites should be included and which sites should be considered
neighbors. As argument to the lattice object, we are required to provide
the interatomic distance.

.. code:: ipython3

    Lat = FCCLattice(IAD=2.770)

Next, we define a **Canvas** object that links Cartesian coordinates to
more abstract graph consisting of sites and neighbors. We incrimentally
construct a Canvas by first introducing a site at the origin of the
coordinate system. Then, we add “two shells” of neighbors, meaning that
we introduce a shell of sites neighboring to the origin (12 for the FCC
lattice) and then introduce another shell of neighbors to that group (42
additional sites, for a total of 55 sites). The lattice object provides
a *getNeighbors* method to identify these neighbors.

.. code:: ipython3

    Canv = Canvas()
    Canv.addLocation(np.array([0,0,0],dtype=float))
    Canv.addShells(2,Lat.getNeighbors)

Finally, we define a list of **Atom** objects to represent the building
blocks of our materials. We then use a **Design** object to represent
the conjunction of a Canvas with a specific arrangement of building
blocks. The Design object will be used to represent the material
decisions made during the solution of material optimization models.

Before applying optimization, we can use the Design object to plot the
sites of the Canvas and ensure that we constructed the intended design
space. We include several parsers to basic crystal structure file
formats such as
`XYZ <https://openbabel.org/docs/dev/FileFormats/XYZ_cartesian_coordinates_format.html>`__,
`PDB <https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html>`__,
`POSCAR <https://cms.mpi.univie.ac.at/vasp/guide/node59.html>`__, and
`CFG <http://li.mit.edu/Archive/Graphics/A/index.html#standard_CFG>`__.

.. code:: ipython3

    Atoms = [Atom('Pt')]
    D = Design(Canv,Atom('Pt'))
    D.toPDB('canvas_sites.pdb')

Building a Model
----------------

To hold the materials information, we create a **MatOptModel** object.
This will hold information about the relevant Canvas, Atoms, and
material conformations that may be present in a system. Additionally, we
define a parameter for the desired size of the cluster which will be
utilized later by several methods.

.. code:: ipython3

    N = 22
    m = MatOptModel(Canv,Atoms)

The MatOptModel additionally hold lists of **MaterialDescriptor**
objects that define the relevant material desriptors. By default,
several universal site descriptors are pre-defined in the model. From
these, all other material descriptors can be defined.

+-------------------------------------+--------------------------------+
| Descriptor                          | Explanation                    |
+=====================================+================================+
| **m.Yik**                           | Presence of a building block   |
|                                     | of type k at site i            |
+-------------------------------------+--------------------------------+
| **m.Yi**                            | Presence of any type of        |
|                                     | building block at site i       |
+-------------------------------------+--------------------------------+
| **m.Xijkl**                         | Presence of a building block   |
|                                     | of type k at site i and a      |
|                                     | building block of type l at    |
|                                     | site j                         |
+-------------------------------------+--------------------------------+
| **m.Xij**                           | Presence of any building block |
|                                     | at site i and any building     |
|                                     | block at site j                |
+-------------------------------------+--------------------------------+
| **m.Cikl**                          | Count of neighbors of type l   |
|                                     | next to a building block of    |
|                                     | type k at site i               |
+-------------------------------------+--------------------------------+
| **m.Ci**                            | Count of any type of neighbors |
|                                     | next to a building block at    |
|                                     | site i                         |
+-------------------------------------+--------------------------------+

User-specified descriptors are defined by **DescriptorRule** objects in
conjunction with **Expr** expression objects. Available expressions
include:

+-----------------------------------------------+-----------------------+
| Expression                                    | Explanation           |
+===============================================+=======================+
| **LinearExpr**                                | Multiplication and    |
|                                               | addition of           |
|                                               | coefficients to       |
|                                               | distinct              |
|                                               | MaterialDescriptors   |
+-----------------------------------------------+-----------------------+
| **SiteCombination**                           | Summation of site     |
|                                               | contributions from    |
|                                               | two sites             |
+-----------------------------------------------+-----------------------+
| **SumNeighborSites**                          | Summation of site     |
|                                               | contributions from    |
|                                               | all neighboring sites |
+-----------------------------------------------+-----------------------+
| **SumNeighborBonds**                          | Summation of bond     |
|                                               | contributions to all  |
|                                               | neighboring sites     |
+-----------------------------------------------+-----------------------+
| **SumSites**                                  | Summation across      |
|                                               | sites                 |
+-----------------------------------------------+-----------------------+
| **SumBonds**                                  | Summation across      |
|                                               | bonds                 |
+-----------------------------------------------+-----------------------+
| **SumSiteTypes**                              | Summation across site |
|                                               | types                 |
+-----------------------------------------------+-----------------------+
| **SumBondTypes**                              | Summation across bond |
|                                               | types                 |
+-----------------------------------------------+-----------------------+
| **SumSitesAndTypes**                          | Summation across      |
|                                               | sites and site types  |
+-----------------------------------------------+-----------------------+
| **SumBondsAndTypes**                          | Summation across      |
|                                               | bonds and bond types  |
+-----------------------------------------------+-----------------------+
| **SumConfs**                                  | Summation across      |
|                                               | conformation types    |
+-----------------------------------------------+-----------------------+
| **SumSitesAndConfs**                          | Summation across      |
|                                               | sites and             |
|                                               | conformation types    |
+-----------------------------------------------+-----------------------+

Several types of DescriptorRules are available.

+-------------------------------------------------+--------------------+
| Rule                                            | Explanation        |
+=================================================+====================+
| **LessThan**                                    | Descriptor less    |
|                                                 | than or equal to   |
|                                                 | an expression      |
+-------------------------------------------------+--------------------+
| **EqualTo**                                     | Descriptor equal   |
|                                                 | to an expression   |
+-------------------------------------------------+--------------------+
| **GreaterThan**                                 | Descriptor greater |
|                                                 | than or equal to   |
|                                                 | an expression      |
+-------------------------------------------------+--------------------+
| **FixedTo**                                     | Descriptor fixed   |
|                                                 | to a scalar value  |
+-------------------------------------------------+--------------------+
| **PiecewiseLinear**                             | Descriptor equal   |
|                                                 | to the evaluation  |
|                                                 | of a piecewise     |
|                                                 | linear function    |
+-------------------------------------------------+--------------------+
| **Implies**                                     | Indicator          |
|                                                 | descriptor that    |
|                                                 | imposes other      |
|                                                 | constraints if     |
|                                                 | equal to 1         |
+-------------------------------------------------+--------------------+
| **NegImplies**                                  | Indicator          |
|                                                 | descriptor that    |
|                                                 | imposes other      |
|                                                 | constraints if     |
|                                                 | equal to 0         |
+-------------------------------------------------+--------------------+
| **ImpliesSiteCombination**                      | Indicator          |
|                                                 | bond-indexed       |
|                                                 | descriptor that    |
|                                                 | imposes            |
|                                                 | constraints on the |
|                                                 | two sites          |
+-------------------------------------------------+--------------------+
| **ImpliesNeighbors**                            | Indicator          |
|                                                 | site-indexed       |
|                                                 | descriptor that    |
|                                                 | imposes            |
|                                                 | constraints on     |
|                                                 | neighboring sites  |
+-------------------------------------------------+--------------------+

From the combination of pre-defined descriptors, expressions, and rules
we can specify a wide variety of other descriptors.

In the context of nanocluster cohesive energy minimization, we would
first like to define the square root of the coordination number. We
achieve this by calling the **addSitesDescriptor** method on
MatOptModel, passing the information necessary to create a
**PiecewiseLinear** rule to correctly define the square root of
coordination at the integer coordination number values. Note that we use
**m.Ci**, the pre-defined basic variable for the count of neighboring
building blocks and equivalent to the coordination number in this
system, as the argument for the piecewise linear function. We use basic
Python lists to express the data for the piecewise linear function
values at integer numbers of coordination.

.. code:: ipython3

    Vals = [sqrt(CN) for CN in range(0,13)]
    BPs = [CN for CN in range(0,13)]
    m.addSitesDescriptor('CNRi',bounds=(0,sqrt(12)),integer=False,
                         rules=PiecewiseLinear(values=Vals,
                                               breakpoints=BPs,
                                               input_desc=m.Ci))

Next, we define a global (i.e., not indexed by sites or bonds)
descriptor for the cohesive energy of the nanocluster. We us a simple
**EqualTo** rule to set this descriptor equal to a normalized sum of
contributions from the square root coordination number descriptor.

.. code:: ipython3

    m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumSites(desc=m.CNRi,
                                                        coefs=(1/(N*sqrt(12))))))

Finally, we create a descriptor for the size of the nanocluster. We set
bounds on this descriptor to effectively constrain the design space to
only include clusters of the desired size, *N*.

.. code:: ipython3

    m.addGlobalDescriptor('Size',bounds=(N,N),
                          rules=EqualTo(SumSites(desc=m.Yi)))

Solving the Model
-----------------

Once all the relevant information in the model is provided, we search
for optimal designs that maximize one of the descriptors. In this
example, we provide the descriptor for coehisver energy as the target
functionality. Additionally, we specify a time limit in seconds as a
keyword argument to the maximize method. For more information, see the
documentation of the maximize function, available in the source code or
by using the Python *help* function.

.. code:: ipython3

    help(MatOptModel.maximize)
    help(MatOptModel.optimize)


.. parsed-literal::

    Help on function maximize in module matopt.opt.mat_modeling:
    
    maximize(self, func, **kwargs)
        Method to maximize a target functionality of the material model.
        
        Args:
        func (MaterialDescriptor/Expr): Material functionality to optimize.
        **kwargs: Arguments to MatOptModel.optimize
        
        Returns:
        (Design/list<Design>) Optimal designs.
        
        See MatOptModel.optimize method for details.
    
    Help on function optimize in module matopt.opt.mat_modeling:
    
    optimize(self, func, sense, nSolns=1, tee=True, disp=1, keepfiles=False, tilim=3600, trelim=None)
        Method to create and optimize the materials design problem.
        
        This method automatically creates a new optimization model every 
        time it is called. Then, it solves the model via Pyomo with the 
        CPLEX solver.
        
        If multiple solutions (called a 'solution pool') are desired, then
        the nSolns argument can be provided and the populate method will 
        be called instead. 
        
        Args:
        func (MaterialDescriptor/Expr): Material functionality to optimize.
        sense (int): flag to indicate the choice to minimize or maximize the
            functionality of interest. 
            Choices: minimize/maximize (Pyomo constants 1,-1 respectively)
        nSolns (int): Optional, number of Design objects to return.
            Default: 1 (See MatOptModel.populate for more information)
        tee (bool): Optional, flag to turn on solver output. 
            Default: True
        disp (int): Optional, flag to control level of MatOpt output.
            Choices: 0: No MatOpt output (other than solver tee)
                     1: MatOpt output for outer level method
                     2: MatOpt output for solution pool & individual solns.
            Default: 1
        keepfiles (bool): Optional, flag to save temporary pyomo files. 
            Default: True
        tilim (float): Optional, solver time limit (in seconds). 
            Default: 3600
        trelim (float): Optional, solver tree memeory limit (in MB).
            Default: None (i.e., Pyomo/CPLEX default)
        
        Returns:
        (Design/list<Design>) Optimal design or designs, depending on the 
            number of solutions requested by argument nSolns.
    


.. code:: ipython3

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
    Logfile '/tmp/tmpyyvw8afg.cplex.log' open.
    CPLEX> New value for absolute mixed integer optimality gap tolerance: 0
    CPLEX> New value for mixed integer optimality gap tolerance: 0
    CPLEX> New value for time limit in seconds: 100
    CPLEX> Problem '/tmp/tmppdm7y76o.pyomo.lp' read.
    Read time = 0.02 sec. (0.18 ticks)
    CPLEX> Problem name         : /tmp/tmppdm7y76o.pyomo.lp
    Objective sense      : Maximize
    Variables            :    1315  [Nneg: 716,  Fix: 1,  Box: 55,  Free: 1,
                                     Binary: 487,  General Integer: 55]
    Objective nonzeros   :       1
    Linear constraints   :    1519  [Less: 1296,  Equal: 223]
      Nonzeros           :    5769
      RHS nonzeros       :     488
    SOS                  :      55  [SOS2: 55, 715 members, all continuous]
    
    Variables            : Min LB: 0.000000         Max UB: 22.00000       
    Objective nonzeros   : Min   : 1.000000         Max   : 1.000000       
    Linear constraints   :
      Nonzeros           : Min   : 0.01312160       Max   : 12.00000       
      RHS nonzeros       : Min   : 1.000000         Max   : 1.000000       
    CPLEX> Tried aggregator 2 times.
    MIP Presolve eliminated 2 rows and 189 columns.
    Aggregator did 55 substitutions.
    Reduced MIP has 1462 rows, 1071 columns, and 5043 nonzeros.
    Reduced MIP has 487 binaries, 0 generals, 55 SOSs, and 0 indicators.
    Presolve time = 0.03 sec. (5.87 ticks)
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
    Root relaxation solution time = 0.03 sec. (40.90 ticks)
    
            Nodes                                         Cuts/
       Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap
    
          0     0        1.2621   271                      1.2621     1179         
    *     0+    0                            0.5285        1.2621           138.81%
    *     0+    0                            0.6055        1.2621           108.42%
          0     0        0.7995   326        0.6055     Cuts: 303     1764   32.04%
    *     0+    0                            0.6310        0.7995            26.71%
          0     0        0.7936   324        0.6310       Cuts: 6     1836   25.77%
          0     0        0.7926   325        0.6310   ZeroHalf: 7     1862   25.61%
          0     0        0.7921   326        0.6310   ZeroHalf: 8     1892   25.53%
          0     0        0.7918   326        0.6310   ZeroHalf: 6     1916   25.49%
    *     0+    0                            0.6639        0.7918            19.26%
          0     0        0.7916   326        0.6639   ZeroHalf: 6     1941   18.90%
          0     0        0.7882   325        0.6639       Cuts: 8     2008   18.71%
          0     0        0.7867   325        0.6639      Cuts: 10     2088   18.49%
    *     0+    0                            0.6790        0.7867            15.87%
    *     0+    0                            0.7087        0.7867            11.01%
          0     0        0.7864   326        0.7087   ZeroHalf: 5     2118   10.97%
          0     0        0.7864   326        0.7087   ZeroHalf: 2     2129   10.96%
          0     0        0.7863   326        0.7087   ZeroHalf: 3     2144   10.95%
    *     0+    0                            0.7251        0.7863             8.44%
          0     2        0.7863   326        0.7251        0.7863     2144    8.44%
    Elapsed time = 0.71 sec. (542.46 ticks, tree = 0.00 MB, solutions = 7)
    *    44+   30                            0.7254        0.7796             7.47%
    *    51    34      integral     0        0.7309        0.7736     6230    5.84%
        165    41        0.7403   174        0.7309        0.7670    13066    4.94%
        781   190        cutoff              0.7309        0.7402    34131    1.27%
    
    Implied bound cuts applied:  55
    Mixed integer rounding cuts applied:  11
    Zero-half cuts applied:  9
    Lift and project cuts applied:  1
    Gomory fractional cuts applied:  7
    
    Root node processing (before b&c):
      Real time             =    0.70 sec. (541.04 ticks)
    Parallel b&c, 8 threads:
      Real time             =    0.61 sec. (639.31 ticks)
      Sync time (average)   =    0.15 sec.
      Wait time (average)   =    0.16 sec.
                              ------------
    Total (root+branch&cut) =    1.32 sec. (1180.35 ticks)
    
    Solution pool: 9 solutions saved.
    
    MIP - Integer optimal solution:  Objective =  7.3092764017e-01
    Solution time =    1.32 sec.  Iterations = 45901  Nodes = 1439
    Deterministic time = 1180.36 ticks  (893.28 ticks/sec)
    
    CPLEX> Incumbent solution written to file '/tmp/tmp480cjg47.cplex.sol'.
    CPLEX> The solver exited normally.
    A feasible and provably optimal solution is available.
    The Design has objective: 0.730927640166201


Processing Results
------------------

If a result is found, we can write it to file and plot with
visualization software. We provide interfaces to several standard
crystal structure file formats, including
`XYZ <https://openbabel.org/docs/dev/FileFormats/XYZ_cartesian_coordinates_format.html>`__,
`PDB <https://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html>`__,
`POSCAR <https://cms.mpi.univie.ac.at/vasp/guide/node59.html>`__, and
`CFG <http://li.mit.edu/Archive/Graphics/A/index.html#standard_CFG>`__.

.. code:: ipython3

    if(D is not None):
        D.toPDB('result.pdb')
