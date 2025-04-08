from pyomo.environ import ConcreteModel, value,units,Constraint
import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Heater
from idaes.core.util.model_statistics import (degrees_of_freedom, unfixed_variables_set)
from idaes.core import MaterialBalanceType
from idaes.core.solvers import get_solver
from HumidAirSurrogate import HAirParameterBlock

vars = ["flow_mass","flow_vol","vol_mass","enth_mass","entr_mass","total_energy_flow"]
indexedVars = ["flow_mol_comp","flow_mass_comp","enth_mol_comp","enth_mass_comp","mass_frac_comp","entr_mol_comp"]
# haha = HAPropsSI("Hha", "T", 416, "P", 5325647, "psi_w", 0.03)
# print(haha)
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = HAirParameterBlock()

# m.fs.heater = Heater(property_package=m.fs.properties, has_pressure_change=False)
# m.fs.heater.inlet.mole_frac_comp[0,"water"].fix(0.01)
# m.fs.heater.inlet.mole_frac_comp[0,"air"].fix(0.99)
# m.fs.heater.inlet.flow_mol[0].fix(1)
# m.fs.heater.inlet.temperature_dry_bulb[0].fix(100+273.15)
# m.fs.heater.inlet.pressure[0].fix(90000)
# m.fs.heater.heat_duty[0].fix(1000)
# # m.fs.heater.deltaP.fix(0)
# flags = m.fs.heater.control_volume.initialize()
# m.fs.heater.control_volume.release_state(flags)
# m.fs.heater.initialize()

solver = get_solver("ipopt")

# print(degrees_of_freedom(m.fs.heater.control_volume.properties_in))
# print(degrees_of_freedom(m.fs.heater.control_volume.properties_out))
# print(degrees_of_freedom(m))
# m.fs.heater.display()

# testSet = unfixed_variables_set(m)
# for comp in testSet:
#     print(comp)
# m.fs.heater.report()


m.fs.sb = m.fs.properties.build_state_block()
m.fs.sb.mole_frac_comp["water"].fix(0.146236327599804)
m.fs.sb.mole_frac_comp["air"].fix(0.853763672400196)
m.fs.sb.flow_mol.fix(1)
m.fs.sb.pressure.fix(269055.991)
m.fs.sb.temperature.fix(499.1001)
m.fs.sb.initialize(outlvl=idaeslog.INFO_HIGH)
print(degrees_of_freedom(m))

solver.solve(m)

# print(value(m.fs.sb.enth_mol))

# print(value(m.fs.sb.enth_mass))

# solver = get_solver("ipopt")
# solver.solve(m)
# # fv.display()
# # print(units.get_units(fv))
# # print(units.get_units(m.fs.sb.enth_mol))
# # print(units.get_units(m.fs.sb.entr_mol))

# print(degrees_of_freedom(m))
m.fs.sb.display()
for var in vars:
#     m.fs.sb.__getattr__(var).display()
#     print(units.get_units(m.fs.sb.__getattr__(var)))
    print(var +": ", value(m.fs.sb.__getattr__(var)))
# for indexedVar in indexedVars:
#     print(indexedVar)
#     for index in m.fs.sb.__getattr__(indexedVar):
#         print(index+": ",value(m.fs.sb.__getattr__(indexedVar)[index]))
# print(value(m.fs.sb.vol_mass))