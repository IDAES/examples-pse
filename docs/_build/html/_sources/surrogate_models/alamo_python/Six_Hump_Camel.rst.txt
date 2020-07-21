.. code:: ipython3

    from idaes.surrogate import alamopy
    import examples
    import pyomo.environ as pyo
    from pyutilib.common import ApplicationError
    import camel6

.. code:: ipython3

    alamo_ran = camel6.main()


.. parsed-literal::

    Updating number of training points from 10 to 10
    Model: {'z1': '  z1 = 3.9999999999999986677324 * x1^2 - 4.0000000000000000000000 * x2^2 - 0.21991890712404245751887E-015 * x1^3 - 2.0999999999999987565502 * x1^4 + 4.0000000000000017763568 * x2^4 + 0.33333333333333303727386 * x1^6 + 0.99999999999999933386619 * x1*x2'}


.. code:: ipython3

    hasBaron = False
    try:
        solver_opt = pyo.SolverFactory('baron')
        hasBaron = solver_opt.available()
    except ApplicationError as err:
        print("Partial test completed. %s"%err)
    
    
    if alamo_ran and hasBaron:
        model = pyo.ConcreteModel()
        opt = pyo.SolverFactory('baron')
        model.x1 = pyo.Var()
        model.x2 = pyo.Var()
        def pyomo_model(model):
            import z1
            return z1.f(model.x1,model.x2)
        model.obj = pyo.Objective(rule = pyomo_model)
        results = opt.solve(model)
        model.solutions.store_to(results)
        print(results)


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
      Termination condition: maxTimeLimit
      Error rc: 0
      Time: 500.21731781959534
    Solution: 
    - number of solutions: 1
      number of solutions displayed: 1
    - Gap: 1e+51
      Status: feasible
      Message: None
      Objective:
        obj:
          Value: -1.031628453489877
      Variable:
        x1:
          Value: 0.08984201313735557
        x2:
          Value: -0.712656403314521
      Constraint: No values
    


