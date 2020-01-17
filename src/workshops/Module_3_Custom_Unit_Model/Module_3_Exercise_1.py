#!/usr/bin/env python
# coding: utf-8

# # Create a Heater

# In this exercise, we will build a heater unit model. The model will assume negligible pressure drop. Most of the code is complete, but you will need to fill in a few key pieces.
# 
# The imports below are complete.

# In[ ]:


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


# Finish the ``make_control_volume`` function by using the control volume to add material, energy, and momentum balances. Try using the documentation for ControlVolume0DBlocks https://idaes-pse.readthedocs.io/en/latest/core/control_volume_0d.html#d-control-volume-class. 
# 
# Hint: For a negligible pressure drop, we can use the total pressure balance defined by ``ControlVolume0DBlock`` without a pressure change.

# In[ ]:


def make_control_volume(unit, name, config):
    if config.dynamic is not False:
        raise ValueError('IdealGasIsentropcCompressor does not support dynamics')
    if config.has_holdup is not False:
        raise ValueError('IdealGasIsentropcCompressor does not support holdup')

    control_volume = ControlVolume0DBlock(default={"property_package": config.property_package,
                                                   "property_package_args": config.property_package_args})

    setattr(unit, name, control_volume)

    control_volume.add_state_blocks(has_phase_equilibrium=config.has_phase_equilibrium)
    ##########
    # Add mass, energy, and momentum balances here
    ##########


# The ``make_config_block`` function is complete.

# In[ ]:


def make_config_block(config):
    config.declare("material_balance_type", ConfigValue(default=MaterialBalanceType.componentPhase, domain=In(MaterialBalanceType)))
    config.declare("energy_balance_type", ConfigValue(default=EnergyBalanceType.enthalpyTotal, domain=In([EnergyBalanceType.enthalpyTotal])))
    config.declare("momentum_balance_type", ConfigValue(default=MomentumBalanceType.pressureTotal, domain=In([MomentumBalanceType.pressureTotal])))
    config.declare("has_phase_equilibrium", ConfigValue(default=False, domain=In([False])))
    config.declare("has_pressure_change", ConfigValue(default=False, domain=In([False])))
    config.declare("property_package", ConfigValue(default=useDefault, domain=is_physical_parameter_block))
    config.declare("property_package_args", ConfigBlock(implicit=True))


# Complete the Heater by defining the ``build`` method.

# In[ ]:


@declare_process_block_class("Heater")
class HeaterData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        # Complete the build method.
        # Don't forget to add a reference to the control 
        # volume 'heat' variable (similar to the compressor 
        # 'work' variable.)
        pass

# The code below is complete. If the Heater is defined correctly, the model below should solve, and the outlet temperature and pressure should be 4.31 (hundred Kelvin) and 0.1 (MPa), respectively.

# In[ ]:


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
opt.options['linear_solver'] = 'mumps'
res = opt.solve(m, tee=False)
print(res.solver.termination_condition)
fs.heater.outlet.display()


# In[ ]:




