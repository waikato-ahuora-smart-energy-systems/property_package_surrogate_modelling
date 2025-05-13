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

#dof = 
print("Dof",degrees_of_freedom(m.fs.state_block))

m.fs.state_block.initialize()

#solve
results = solver.solve(m, tee=tee)


# Access and initialize each property
print("Flow Mass:", value(m.fs.state_block.flow_mass))
print("Flow Mass Comp:", {k: value(v) for k, v in m.fs.state_block.flow_mass_comp.items()})
print("Flow Vol:", value(m.fs.state_block.flow_vol))
print("Enthalpy Molar:", value(m.fs.state_block.enth_mol))
print("Enthalpy Mass:", value(m.fs.state_block.enth_mass))
print("Mole Fraction Comp:", {k: value(v) for k, v in m.fs.state_block.mole_frac_comp.items()})
print("Mass Fraction Comp:", {k: value(v) for k, v in m.fs.state_block.mass_frac_comp.items()})
print("Entropy Molar:", value(m.fs.state_block.entr_mol))
print("Entropy Mass:", value(m.fs.state_block.entr_mass))
print("Flow Mass Phase:", {k: value(v) for k, v in m.fs.state_block.flow_mass_phase.items()})
print("Flow Molar Phase:", {k: value(v) for k, v in m.fs.state_block.flow_mol_phase.items()})
print("Mole Fraction Phase Comp:", {k: value(v) for k, v in m.fs.state_block.mole_frac_phase_comp.items()})
print("Mass Fraction Phase Comp:", {k: value(v) for k, v in m.fs.state_block.mass_frac_phase_comp.items()})
print("Phase Fraction Vap:", value(m.fs.state_block.phase_frac_vap))
print("Phase Fraction Liq:", value(m.fs.state_block.phase_frac_liq))
print("Spec Vol Mol Phase Liq:", value(m.fs.state_block.spec_vol_mol_phase["Liq"]))
print("Spec Vol Mol Phase Vap:", value(m.fs.state_block.spec_vol_mol_phase["Vap"]))
print("Enth molar Phase Vap:", value(m.fs.state_block.enth_mol_vap))
print("Enth Molar Phase Liq:", value(m.fs.state_block.enth_mol_phase["Liq"]))
print("Enth Molar Phase Vap:", value(m.fs.state_block.enth_mol_phase["Vap"]))
print("Enth Mass Phase Liq:", value(m.fs.state_block.enth_mass_phase["Liq"]))
print("Enth Mass Phase Vap:", value(m.fs.state_block.enth_mass_phase["Vap"]))
print("Entr Molar Phase Liq:", value(m.fs.state_block.entr_mol_phase["Liq"]))
print("Entr Molar Phase Vap:", value(m.fs.state_block.entr_mol_phase["Vap"]))
print("Entr Mass Phase liq water:", value(m.fs.state_block.entr_mass_phase_comp["Liq", "water"]))
print("Entr Mass Phase vap water:", value(m.fs.state_block.entr_mass_phase_comp["Vap", "water"]))
print("Entr Mass Phase liq air:", value(m.fs.state_block.entr_mass_phase_comp["Liq", "air"]))
print("Entr Mass Phase vap air:", value(m.fs.state_block.entr_mass_phase_comp["Vap", "air"]))
print("Enth Mass Phase vap air:", value(m.fs.state_block.enth_mass_phase_comp["Vap", "air"]))
print("Enth Mass Phase liq air:", value(m.fs.state_block.enth_mass_phase_comp["Liq", "air"]))
print("Enth Mass Phase vap water:", value(m.fs.state_block.enth_mass_phase_comp["Vap", "water"]))
print("Enth Mass Phase liq water:", value(m.fs.state_block.enth_mass_phase_comp["Liq", "water"]))
print("Flow Mass Phase vap air:", value(m.fs.state_block.flow_mass_phase_comp["Vap", "air"]))
m.fs.state_block.flow_mass
m.fs.state_block.flow_mass_comp
print(value(m.fs.state_block.flow_mol))
print(value(m.fs.state_block.phase_frac_vap))
print(value(m.fs.state_block.spec_vol_mol_phase["Liq"]))
print(value(m.fs.state_block.phase_frac_liq))
print(value(m.fs.state_block.vol_mol_vap))
print(value(m.fs.state_block.flow_vol))
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




m.fs.state_block.display()