
.. code:: ipython2

    from idaes.surrogate import alamopy

.. code:: ipython2

    execfile('camel6.py')


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-2-4f88f14d452a> in <module>
    ----> 1 execfile('camel6.py')
    

    NameError: name 'execfile' is not defined


.. code:: ipython2

    import pyomo.environ as pyo
    model = pyo.ConcreteModel()
    opt = pyo.SolverFactory('baron')
    model.x1 = pyo.Var()
    model.x2 = pyo.Var()
    def pyomo_model(model):
        import cam6alm
        return cam6alm.f(model.x1,model.x2)
    model.obj = pyo.Objective(rule = pyomo_model)
    results = opt.solve(model)

.. code:: ipython2

    model.solutions.store_to(results)
    print results


.. parsed-literal::

    
    Problem: 
    - Name: problem
      Lower bound: -1e+51
      Upper bound: -1.03162845349
      Number of objectives: 1
      Number of constraints: 1
      Number of variables: 3
      Sense: unknown
    Solver: 
    - Status: ok
      Termination condition: unknown
      Error rc: 0
      Time: 0.753000020981
    Solution: 
    - number of solutions: 1
      number of solutions displayed: 1
    - Gap: 1e+51
      Status: feasible
      Message: None
      Objective:
        obj:
          Value: -1.03162845349
      Variable:
        x1:
          Value: 0.0898426337843
        x2:
          Value: -0.712656468875
      Constraint: No values
    


