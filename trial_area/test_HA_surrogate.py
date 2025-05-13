from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    Param,
    Constraint,
    SolverFactory,
    TransformationFactory,
    TerminationCondition,
    value,
    Expression,
    minimize,
    units,
)
from pyomo.network import Arc, SequentialDecomposition

# Import IDAES libraries
from idaes.core import FlowsheetBlock, UnitModelBlockData
from idaes.models.unit_models import (
    Mixer,
    MomentumMixingType,
    PressureChanger,
    Heater,
    Separator,
    HeatExchanger,
)
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from HumidAirSurrogate import HAirParameterBlock

import idaes.logger as idaeslog


# Setup solver and options
solver = SolverFactory("ipopt")
outlvl = 0
tee = True

# Set up concrete model
m = ConcreteModel()

# Create a flowsheet block
m.fs = FlowsheetBlock(dynamic=False)

# Create the properties param block
m.fs.properties = HAirParameterBlock()

# Create a state block for the humid air properties
m.fs.state_block = m.fs.properties.build_state_block(
    defined_state= True
)

# Example variables for the state block
m.fs.state_block.temperature.fix(300)
m.fs.state_block.pressure.fix(101325)
m.fs.state_block.flow_mol.fix(1)
m.fs.state_block.mole_frac_comp["water"].fix(0.01)
m.fs.state_block.mole_frac_comp["air"].fix(0.99)




# Access and initialize each property
m.fs.state_block.flow_mass
m.fs.state_block.flow_mass_comp
m.fs.state_block.flow_vol
m.fs.state_block.enth_mol
m.fs.state_block.enth_mass
m.fs.state_block.mole_frac_comp
m.fs.state_block.mass_frac_comp
m.fs.state_block.entr_mol
m.fs.state_block.entr_mass
m.fs.state_block.flow_mass_phase
m.fs.state_block.flow_mol_phase
m.fs.state_block.mole_frac_phase_comp
m.fs.state_block.mass_frac_phase_comp

#dof = 
print("Dof",degrees_of_freedom(m.fs.state_block))

m.fs.state_block.initialize()

#solve
results = solver.solve(m, tee=tee)


m.fs.state_block.display()