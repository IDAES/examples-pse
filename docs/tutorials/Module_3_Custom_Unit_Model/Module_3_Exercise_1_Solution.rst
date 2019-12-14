Create a Heater
===============

.. code:: ipython3

    import pyomo.environ as pe
    from pyomo.common.config import ConfigBlock, ConfigValue, In
    from idaes.core import (ControlVolume0DBlock,
                            declare_process_block_class,
                            EnergyBalanceType,
                            MomentumBalanceType,
                            MaterialBalanceType,
                            UnitModelBlockData,
                            useDefault,
                            FlowsheetBlock)
    from idaes.core.util.config import is_physical_parameter_block
    from methanol_param_VLE import PhysicalParameterBlock
    from idaes.core.util.misc import add_object_reference

.. code:: ipython3

    def make_control_volume(unit, name, config):
        if config.dynamic is not False:
            raise ValueError('IdealGasIsentropcCompressor does not support dynamics')
        if config.has_holdup is not False:
            raise ValueError('IdealGasIsentropcCompressor does not support holdup')
    
        control_volume = ControlVolume0DBlock(default={"property_package": config.property_package,
                                                       "property_package_args": config.property_package_args})
    
        setattr(unit, name, control_volume)
    
        control_volume.add_state_blocks(has_phase_equilibrium=config.has_phase_equilibrium)
        control_volume.add_material_balances(balance_type=config.material_balance_type,
                                             has_phase_equilibrium=config.has_phase_equilibrium)
        control_volume.add_total_enthalpy_balances(has_heat_of_reaction=False, 
                                                   has_heat_transfer=True, 
                                                   has_work_transfer=False)
        control_volume.add_total_pressure_balances(has_pressure_change=False)

.. code:: ipython3

    def make_config_block(config):
        config.declare("material_balance_type",
            ConfigValue(default=MaterialBalanceType.componentPhase, domain=In(MaterialBalanceType)))
        config.declare("energy_balance_type",
            ConfigValue(default=EnergyBalanceType.enthalpyTotal, domain=In([EnergyBalanceType.enthalpyTotal])))
        config.declare("momentum_balance_type",
            ConfigValue(default=MomentumBalanceType.pressureTotal, domain=In([MomentumBalanceType.pressureTotal])))
        config.declare("has_phase_equilibrium",
            ConfigValue(default=False, domain=In([False])))
        config.declare("has_pressure_change",
            ConfigValue(default=False, domain=In([False])))
        config.declare("property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block))
        config.declare("property_package_args",
            ConfigBlock(implicit=True))

.. code:: ipython3

    @declare_process_block_class("Heater")
    class HeaterData(UnitModelBlockData):
        CONFIG = UnitModelBlockData.CONFIG()
        make_config_block(CONFIG)
    
        def build(self):
            super(HeaterData, self).build()
    
            make_control_volume(self, "control_volume", self.config)
    
            self.add_inlet_port()
            self.add_outlet_port()
            
            add_object_reference(self, 'heat', self.control_volume.heat[0.0])

.. code:: ipython3

    m = pe.ConcreteModel()
    m.fs = fs = FlowsheetBlock(default={"dynamic": False})
    fs.properties = props = PhysicalParameterBlock(default={'Cp': 0.038056, 'valid_phase': 'Vap'})
    
    fs.heater = Heater(default={"property_package": props, 'has_phase_equilibrium': False})
    fs.heater.inlet.flow_mol.fix(1)
    fs.heater.inlet.mole_frac_comp[0, 'CH3OH'].fix(0.25)
    fs.heater.inlet.mole_frac_comp[0, 'CH4'].fix(0.25)
    fs.heater.inlet.mole_frac_comp[0, 'H2'].fix(0.25)
    fs.heater.inlet.mole_frac_comp[0, 'CO'].fix(0.25)
    fs.heater.inlet.pressure.fix(0.1)
    fs.heater.inlet.temperature.fix(3)
    fs.heater.heat.fix(5)
    
    opt = pe.SolverFactory('ipopt')
    res = opt.solve(m, tee=True)
    print(res.solver.termination_condition)
    fs.heater.outlet.display()


.. parsed-literal::

    Ipopt 3.12.13: 
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    
    This version of Ipopt was compiled from source code available at
        https://github.com/IDAES/Ipopt as part of the Institute for the Design of
        Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
        Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    
    This version of Ipopt was compiled using HSL, a collection of Fortran codes
        for large-scale scientific computation.  All technical papers, sales and
        publicity material resulting from use of the HSL codes within IPOPT must
        contain the following acknowledgement:
            HSL, a collection of Fortran codes for large-scale scientific
            computation. See http://www.hsl.rl.ac.uk.
    ******************************************************************************
    
    This is Ipopt version 3.12.13, running with linear solver ma27.
    
    Number of nonzeros in equality constraint Jacobian...:       51
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       13
    
    Total number of variables............................:       17
                         variables with only lower bounds:        5
                    variables with lower and upper bounds:       12
                         variables with only upper bounds:        0
    Total number of equality constraints.................:       17
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  0.0000000e+00 5.00e-01 1.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 5.00e-06 5.26e+00  -1.0 2.63e+00    -  5.37e-01 1.00e+00h  1
       2  0.0000000e+00 0.00e+00 2.22e+00  -1.7 1.31e+00    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 2
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   0.0000000000000000e+00    0.0000000000000000e+00
    
    
    Number of objective function evaluations             = 3
    Number of objective gradient evaluations             = 3
    Number of equality constraint evaluations            = 3
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 3
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 2
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.000
    Total CPU secs in NLP function evaluations           =      0.000
    
    EXIT: Optimal Solution Found.
    optimal
    outlet : Size=1
        Key  : Name           : Value
        None :       flow_mol : {0.0: 1.0}
             : mole_frac_comp : {(0.0, 'CH3OH'): 0.25, (0.0, 'CH4'): 0.25, (0.0, 'CO'): 0.25, (0.0, 'H2'): 0.25}
             :       pressure : {0.0: 0.1}
             :    temperature : {0.0: 4.313853268866933}


.. code:: ipython3

    # For testing purposes
    from pyomo.environ import TerminationCondition
    assert res.solver.termination_condition == TerminationCondition.optimal


