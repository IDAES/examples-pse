.. code:: ipython3

    from idaes.surrogate import alamopy
    import examples
    import pyomo.environ as pyo
    import camel6

.. code:: ipython3

    alamo_ran = camel6.main()

.. code:: ipython3

    if alamo_ran:
        model = pyo.ConcreteModel()
        opt = pyo.SolverFactory('baron')
        model.x1 = pyo.Var()
        model.x2 = pyo.Var()
        def pyomo_model(model):
            import cam6alm
            return cam6alm.f(model.x1,model.x2)
        model.obj = pyo.Objective(rule = pyomo_model)
        results = opt.solve(model)
        model.solutions.store_to(results)
        print(results)

