from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value, \
    Expression, Objective, minimize
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock

import pyomo.contrib.parmest.parmest as parmest

# Set tag level to see output logs; here we want to see log messages
# down to the properties level
idaeslog.add_log_tag("properties")

# Create model function
def NRTL_model(data):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                                 ('Liq', 'Vap'),
                                                 "activity_coeff_model":
                                                 'NRTL'})
    m.fs.state_block = m.fs.properties.state_block_class(
        default={"parameters": m.fs.properties,
                 "defined_state": True})

    # Initialize at a certain inlet condition
    m.fs.state_block.flow_mol.fix(1)
    m.fs.state_block.temperature.fix(368)
    m.fs.state_block.pressure.fix(101325)
    m.fs.state_block.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block.mole_frac_comp["toluene"].fix(0.5)

    # Fix NRTL specific variables
    # alpha values (set at 0.3)
    m.fs.state_block.\
        alpha["benzene", "benzene"].fix(0)
    m.fs.state_block.\
        alpha["benzene", "toluene"].fix(0.3)
    m.fs.state_block.\
        alpha["toluene", "toluene"].fix(0)
    m.fs.state_block.\
        alpha["toluene", "benzene"].fix(0.3)

    # initial tau values
    m.fs.state_block.\
        tau["benzene", "benzene"].fix(0)
    m.fs.state_block.\
        tau["benzene", "toluene"].fix(0.1690)
    m.fs.state_block.\
        tau["toluene", "toluene"].fix(0)
    m.fs.state_block.\
        tau["toluene", "benzene"].fix(-0.1559)

    # Initialize the flash unit
    m.fs.state_block.initialize(outlvl=idaeslog.INFO)

    # Fix at actual temperature
    m.fs.state_block.temperature.fix(float(data["temperature"]))

    # Set bounds on variables to be estimated
    m.fs.state_block.\
        tau["benzene", "toluene"].setlb(-5)
    m.fs.state_block.\
        tau["benzene", "toluene"].setub(5)

    m.fs.state_block.\
        tau["toluene", "benzene"].setlb(-5)
    m.fs.state_block.\
        tau["toluene", "benzene"].setub(5)

    # Return initialized flash model
    return m

# Parameter Estimation using parmest

# Vars to estimate
variable_name = ["fs.state_block.tau['benzene', 'toluene']",
                 "fs.state_block.tau['toluene', 'benzene']"]

# List of dictionaries
data = [{"temperature": 367, "vap_benzene": 0.662038, "liq_benzene": 0.441008},
        {"temperature": 368, "vap_benzene": 0.631425, "liq_benzene": 0.409061},
        {"temperature": 369, "vap_benzene": 0.599865, "liq_benzene": 0.378040},
        {"temperature": 370, "vap_benzene": 0.567323, "liq_benzene": 0.347899},
        {"temperature": 371, "vap_benzene": 0.533776, "liq_benzene": 0.318596}]

# Create expression to compute the sum of squared error
def SSE(m, data):
    expr = ((float(data["vap_benzene"]) -
             m.fs.state_block.mole_frac_phase_comp["Vap", "benzene"])**2 +
            (float(data["liq_benzene"]) -
             m.fs.state_block.mole_frac_phase_comp["Liq", "benzene"])**2)
    return expr

pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True)

SSE, parameters = pest.theta_est()
print(SSE)
print(parameters)
