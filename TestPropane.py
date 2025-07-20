import matplotlib.pyplot as plt
#Importing required pyomo and idaes components
from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    Objective,
    SolverFactory,
    TransformationFactory,
    value,
    units,
)
from pyomo.network import Arc, SequentialDecomposition

from pyomo.util.infeasible import (
    log_infeasible_constraints,
    log_infeasible_bounds)
import logging

from idaes.core import FlowsheetBlock,MaterialBalanceType

from idaes.models.unit_models import (
    PressureChanger,
    Mixer,
    Separator as Splitter,
    Heater,
    Compressor,
    HeatExchanger
)
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.models.unit_models.heat_exchanger import HX0DInitializer
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
from idaes.core.util.model_statistics import degrees_of_freedom

# Import idaes logger to set output levels
import idaes.logger as idaeslog

from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    HelmholtzParameterBlockData,
    AmountBasis
)

from idaes.models.properties.modular_properties.coolprop.coolprop_wrapper import (
    CoolPropWrapper,
    CoolPropExpressionError,
    CoolPropPropertyError,
)
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
#Constructing the Flowsheet

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
logging.basicConfig(filename='PyomoLog.log', encoding='utf-8', level=logging.INFO)

#m.fs.propertiesR = GenericParameterBlock(**configuration)

m.fs.propertiesRC = HelmholtzParameterBlock(
  pure_component="Propane",
  phase_presentation=PhaseType.MIX,
  state_vars=StateVars.PH,
)
#m.fs.propertiesRC.default_enthalpy_mol_bounds = (-1000, 350000)

# Set up the flowsheet

m.fs.Exp = PressureChanger(
    property_package=m.fs.propertiesRC,
    compressor=False,
    thermodynamic_assumption=ThermodynamicAssumption.adiabatic,
)

m.fs.Exp.inlet.flow_mol[0].fix(100)
m.fs.Exp.inlet.enth_mol[0].fix(m.fs.propertiesRC.htpx(p=2000000 * units.Pa, T=298.15*units.K))
#m.fs.Exp.inlet.temperature[0].fix(298.15)
m.fs.Exp.inlet.pressure[0].fix(2000000)
#m.fs.Exp.inlet.vapor_frac[0].fix(0)

m.fs.Exp.outlet.pressure[0].fix(500000)

#Degrees of freedom check
print(degrees_of_freedom(m))

# Initialize the model
m.fs.Exp.initialize()
m.fs.Exp.report()

# Solve the model
solver = SolverFactory("ipopt")
results = solver.solve(m, tee=True)


m.fs.Exp.outlet.pressure.unfix()
@m.fs.Exp.control_volume.properties_out[0].Constraint()
def temp_constraint(b):
    return b.temperature == 1.74 + 273.15
      
#m.fs.Exp.outlet.temperature.fix(1.74+273.15)

results = solver.solve(m, tee=True)
m.fs.Exp.report()

m.fs.Exp.control_volume.properties_out[0].temp_constraint.deactivate()
m.fs.Exp.initialize()
m.fs.Exp.report()

m.fs.Exp.control_volume.properties_out[0].temp_constraint.activate()
results = solver.solve(m, tee=True)
# m.fs.Exp.control_volume.properties_out[0].display()
m.fs.Exp.report()
m.fs.Exp.control_volume.properties_out[0].pressure.set_value(50000)
results = solver.solve(m, tee=True)

m.fs.Exp.report()
print(degrees_of_freedom(m))