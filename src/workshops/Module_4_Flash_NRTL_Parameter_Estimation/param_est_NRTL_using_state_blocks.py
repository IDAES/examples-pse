from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock

import pyomo.contrib.parmest.parmest as parmest
import pandas as pd
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
    m.fs.properties.\
        alpha["benzene", "benzene"].fix(0)
    m.fs.properties.\
        alpha["benzene", "toluene"].fix(0.3)
    m.fs.properties.\
        alpha["toluene", "toluene"].fix(0)
    m.fs.properties.\
        alpha["toluene", "benzene"].fix(0.3)

    # initial tau values
    m.fs.properties.\
        tau["benzene", "benzene"].fix(0)
    m.fs.properties.\
        tau["benzene", "toluene"].fix(0.1690)
    m.fs.properties.\
        tau["toluene", "toluene"].fix(0)
    m.fs.properties.\
        tau["toluene", "benzene"].fix(-0.1559)

    # Initialize the flash unit
    m.fs.state_block.initialize(outlvl=idaeslog.INFO)

    # Fix at actual temperature
    m.fs.state_block.temperature.fix(float(data["temperature"]))

    # Set bounds on variables to be estimated
    m.fs.properties.\
        tau["benzene", "toluene"].setlb(-5)
    m.fs.properties.\
        tau["benzene", "toluene"].setub(5)

    m.fs.properties.\
        tau["toluene", "benzene"].setlb(-5)
    m.fs.properties.\
        tau["toluene", "benzene"].setub(5)

    # Return initialized flash model
    return m

# Parameter Estimation using parmest

# Vars to estimate
variable_name = ["fs.properties.tau['benzene', 'toluene']",
                 "fs.properties.tau['toluene', 'benzene']"]

# Load data
data = pd.read_csv('BT_NRTL_dataset.csv')

# Create expression to compute the sum of squared error
def SSE(m, data):
    expr = ((float(data["vap_benzene"]) -
             m.fs.state_block.mole_frac_phase_comp["Vap", "benzene"])**2 +
            (float(data["liq_benzene"]) -
             m.fs.state_block.mole_frac_phase_comp["Liq", "benzene"])**2)
    return expr*1E4

# Initialize a parameter estimation object
pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True)

# Run parameter estimation using all data
obj_value, parameters = pest.theta_est()
print(obj_value)
print(parameters)

# Run parameter estimation using bootstrap resample of the data (8 samples), 
# plot results along with confidence regions
bootstrap_theta = pest.theta_est_bootstrap(8)
print(bootstrap_theta)
parmest.pairwise_plot(bootstrap_theta, alpha=0.75, distributions=['Rect', 'MVN'])

# Run parameter estimation using leave-N-out samples (leave 2 out, 8 samples), 
# plot results along with confidence regions
lNo_theta = pest.theta_est_leaveNout(2, 8)
print(lNo_theta)
parmest.pairwise_plot(lNo_theta, alpha=0.75, distributions=['Rect', 'MVN'])

# Run a confidence region test and plot data points that 
# fall within an alpha confidence region
test_results = pest.confidence_region_test(bootstrap_theta, distribution='MVN', 
                                           alphas=[0.75, 0.9])
print(test_results)
parmest.pairwise_plot(test_results, alpha=0.75)

