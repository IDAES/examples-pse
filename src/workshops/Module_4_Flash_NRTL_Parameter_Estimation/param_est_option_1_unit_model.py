from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value, \
    Expression, Objective, minimize
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock


# Set tag level to see output logs; here we want to see log messages down to the properties level
idaeslog.add_log_tag("properties")

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             ('Liq', 'Vap'),
                                             "activity_coeff_model":
                                             'NRTL'})

from idaes.generic_models.unit_models import Flash


m.fs.flash = Flash(default={"property_package": m.fs.properties})

from idaes.core.util.model_statistics import degrees_of_freedom

# Todo: print the degrees of freedom for your model
print("Degrees of Freedom =", degrees_of_freedom(m))


# Inlet specifications given above
m.fs.flash.inlet.flow_mol.fix(1)
m.fs.flash.inlet.temperature.fix(368)
m.fs.flash.inlet.pressure.fix(101325)
m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

# Todo: add code for the 2 flash unit specifications given above
m.fs.flash.heat_duty.fix(0)
m.fs.flash.deltaP.fix(0)

# Fix NRTL specific
# alpha values (set at 0.3)
m.fs.flash.control_volume.properties_in[0].alpha["benzene", "benzene"].fix(0)
m.fs.flash.control_volume.properties_in[0].alpha["benzene", "toluene"].fix(0.3)
m.fs.flash.control_volume.properties_in[0].alpha["toluene", "toluene"].fix(0)
m.fs.flash.control_volume.properties_in[0].alpha["toluene", "benzene"].fix(0.3)

m.fs.flash.control_volume.properties_out[0].\
    alpha["benzene", "benzene"].fix(0)
m.fs.flash.control_volume.properties_out[0].\
    alpha["benzene", "toluene"].fix(0.3)
m.fs.flash.control_volume.properties_out[0].\
    alpha["toluene", "toluene"].fix(0)
m.fs.flash.control_volume.properties_out[0].\
    alpha["toluene", "benzene"].fix(0.3)

# tau values
m.fs.flash.control_volume.properties_in[0].\
    tau["benzene", "benzene"].fix(0)
m.fs.flash.control_volume.properties_in[0].\
    tau["benzene", "toluene"].fix(0.3)
m.fs.flash.control_volume.properties_in[0].\
    tau["toluene", "toluene"].fix(0)
m.fs.flash.control_volume.properties_in[0].\
    tau["toluene", "benzene"].fix(-0.1)

m.fs.flash.control_volume.properties_out[0].\
    tau["benzene", "benzene"].fix(0)
m.fs.flash.control_volume.properties_out[0].\
    tau["benzene", "toluene"].fix(0.2)
m.fs.flash.control_volume.properties_out[0].\
    tau["toluene", "toluene"].fix(0)
m.fs.flash.control_volume.properties_out[0].\
    tau["toluene", "benzene"].fix(-0.1)

# Todo: print the degrees of freedom for your model
print("Degrees of Freedom =", degrees_of_freedom(m))

# Todo: initialize the flash unit
m.fs.flash.initialize(outlvl=idaeslog.INFO)

# Todo: create the ipopt solver
solver = SolverFactory('ipopt')

# Todo: solve the model
status = solver.solve(m, tee=True)

m.fs.flash.report()
m.fs.flash.control_volume.properties_out[0].temperature_bubble.display()

# Parameter Estimation problem setup
VLE_data = {"vap_benzene": 0.631425, "liq_benzene": 0.40906}
m.fs.sse = Expression(expr=(VLE_data["vap_benzene"] -
                            m.fs.flash.vap_outlet.
                            mole_frac_comp[0, "benzene"])**2 +
                      (VLE_data["liq_benzene"] -
                       m.fs.flash.liq_outlet.
                       mole_frac_comp[0, "benzene"])**2)

m.fs.param_obj = Objective(expr=m.fs.sse, sense=minimize)

m.fs.flash.control_volume.properties_in[0].tau["benzene", "toluene"].unfix()
m.fs.flash.control_volume.properties_in[0].tau["toluene", "benzene"].unfix()

m.fs.flash.control_volume.properties_out[0].tau["benzene", "toluene"].unfix()
m.fs.flash.control_volume.properties_out[0].tau["toluene", "benzene"].unfix()

# Add constraints equating tau from properties_in and from properties_out
m.fs.c1 = Constraint(expr=m.fs.flash.control_volume.
                     properties_in[0].tau["benzene", "toluene"] ==
                     m.fs.flash.control_volume.
                     properties_out[0].tau["benzene", "toluene"])

m.fs.c2 = Constraint(expr=m.fs.flash.control_volume.
                     properties_in[0].tau["toluene", "benzene"] ==
                     m.fs.flash.control_volume.
                     properties_out[0].tau["toluene", "benzene"])

m.fs.flash.control_volume.properties_in[0].\
    tau["benzene", "toluene"].setlb(0)
m.fs.flash.control_volume.properties_in[0].\
    tau["benzene", "toluene"].setub(1)

m.fs.flash.control_volume.properties_in[0].\
    tau["toluene", "benzene"].setlb(-1)
m.fs.flash.control_volume.properties_in[0].\
    tau["toluene", "benzene"].setub(1)


print(degrees_of_freedom(m))

status = solver.solve(m, tee=True)

m.fs.flash.report()
m.fs.flash.control_volume.properties_in[0].tau.display()
