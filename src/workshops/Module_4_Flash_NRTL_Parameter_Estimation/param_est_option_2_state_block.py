from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value, \
    Expression, Objective, minimize
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Set tag level to see output logs; here we want to see log messages
# down to the properties level
idaeslog.add_log_tag("properties")

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             ('Liq', 'Vap'),
                                             "activity_coeff_model":
                                             'NRTL'})
m.fs.state_block = m.fs.properties.state_block_class(
    default={"parameters": m.fs.properties,
             "defined_state": True})


# Todo: print the degrees of freedom for your model
print("Degrees of Freedom =", degrees_of_freedom(m))

m.fs.state_block.flow_mol.fix(1)
m.fs.state_block.temperature.fix(368)
m.fs.state_block.pressure.fix(101325)
m.fs.state_block.mole_frac_comp["benzene"].fix(0.5)
m.fs.state_block.mole_frac_comp["toluene"].fix(0.5)

print("Degrees of Freedom =", degrees_of_freedom(m))


# Fix NRTL specific
# alpha values (set at 0.3)
m.fs.properties.alpha["benzene", "benzene"].fix(0)
m.fs.properties.alpha["benzene", "toluene"].fix(0.3)
m.fs.properties.alpha["toluene", "toluene"].fix(0)
m.fs.properties.alpha["toluene", "benzene"].fix(0.3)

# tau values
m.fs.properties.tau["benzene", "benzene"].fix(0)
m.fs.properties.tau["benzene", "toluene"].fix(0.1690)
m.fs.properties.tau["toluene", "toluene"].fix(0)
m.fs.properties.tau["toluene", "benzene"].fix(-0.1559)

# Todo: print the degrees of freedom for your model
print("Degrees of Freedom =", degrees_of_freedom(m))

# Todo: initialize the flash unit
m.fs.state_block.initialize(outlvl=idaeslog.INFO)

# Parameter Estimation
VLE_data = {"vap_benzene": 0.631425, "liq_benzene": 0.40906}

# Create expression to compute the sum of squared errors
m.fs.sse = Expression(expr=(VLE_data["vap_benzene"] -
                            m.fs.state_block.
                            mole_frac_phase_comp["Vap", "benzene"])**2 +
                      (VLE_data["liq_benzene"] -
                       m.fs.state_block.
                       mole_frac_phase_comp["Liq", "benzene"])**2)

m.fs.param_obj = Objective(expr=m.fs.sse, sense=minimize)

m.fs.properties.tau["benzene", "toluene"].unfix()
m.fs.properties.tau["toluene", "benzene"].unfix()

m.fs.properties.tau["benzene", "toluene"].setlb(-5)
m.fs.properties.tau["benzene", "toluene"].setub(5)

m.fs.properties.tau["toluene", "benzene"].setlb(-5)
m.fs.properties.tau["toluene", "benzene"].setub(5)


print(degrees_of_freedom(m))

solver = SolverFactory("ipopt")
status = solver.solve(m, tee=True)


m.fs.properties.tau.display()

m.fs.state_block.mole_frac_phase_comp.display()
