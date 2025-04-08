import pytest
from pyomo.environ import ConcreteModel, value, SolverFactory
from idaes.core import FlowsheetBlock
from humid_air_stateblock import HAirParameterBlock

@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = HAirParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)

    s = m.fs.state[0]
    s.flow_mol.fix(1.0)
    s.temperature.fix(300)
    s.pressure.fix(101325)
    s.mole_frac_comp['water'].fix(0.02)
    s.mole_frac_comp['air'].fix(0.98)

    return m

def test_state_vars_exist(model):
    s = model.fs.state[0]
    for attr in ["flow_mol", "temperature", "pressure", "mole_frac_comp"]:
        assert hasattr(s, attr)

# def test_surrogate_outputs_exist(model):
#     s = model.fs.state[0]
#     props = [
#         "relative_humidity", "temperature_wet_bulb",
#         "enth_mol_gas_phase", "entr_mol_gas_phase", "vol_mol_gas_phase",
#         "enth_mol_liq_phase", "entr_mol_liq_phase", "vol_mol_liq_phase",
#         "mol_frac_gas_phase", "mole_frac_1_gas_phase", "mole_frac_1_liq_phase",
#         "enth_mol", "entr_mol", "vol_mol"
#     ]
#     for prop in props:
#         assert hasattr(s, prop)

# def test_model_solves(model):
#     solver = SolverFactory("ipopt")
#     result = solver.solve(model, tee=False)
#     assert result.solver.termination_condition.name == "optimal"
#     assert result.solver.status == result.solver.status.ok

# def test_outputs_are_finite(model):
#     s = model.fs.state[0]
#     value_props = [
#         s.relative_humidity,
#         s.temperature_wet_bulb,
#         s.enth_mol,
#         s.entr_mol,
#         s.vol_mol,
#     ]
#     for prop in value_props:
#         assert value(prop) != pytest.approx(float("inf"))
#         assert value(prop) != pytest.approx(float("nan"))

# def test_bounds(model):
#     s = model.fs.state[0]
#     assert 0 <= value(s.relative_humidity) <= 1
#     assert 200 < value(s.temperature_wet_bulb) < value(s.temperature)
#     assert value(s.vol_mol) > 0
#     assert value(s.enth_mol) > 0
#     assert value(s.entr_mol) > 0
