from pyomo.environ import units, Var, Expression, Reals, NonNegativeReals
from idaes.core import (
    declare_process_block_class, PhysicalParameterBlock,
    StateBlockData, StateBlock, VaporPhase, LiquidPhase,
    Component, MaterialBalanceType, EnergyBalanceType, PhaseType
)
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
import os
import warnings

@declare_process_block_class("HAirParameterBlock")
class HAirParameterBlockData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.Vap = VaporPhase()
        self.Liq = LiquidPhase()

        self.air = Component()
        self.water = Component()

        self.mw_comp = self.define_metadata().add_default_units({
            "mass": units.kg,
            "amount": units.mol
        })

        self.mw_comp = self.define_parameter(
            self.component_list,
            initialize={"water": 0.01801528, "air": 0.02896},
            units=units.kg / units.mol
        )

        self._state_block_class = HAirStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            "flow_mol": {"method": None, "units": units.mol / units.s},
            "temperature": {"method": None, "units": units.K},
            "pressure": {"method": None, "units": units.Pa},
            "mole_frac_comp": {"method": None},
            "enth_mol_phase": {"method": None, "units": units.J / units.mol},
            "entr_mol_phase": {"method": None, "units": units.J / units.mol / units.K},
            "vol_mol_phase": {"method": None, "units": units.m**3 / units.mol}
        })

        obj.define_custom_properties({
            "relative_humidity": {"method": None},
            "temperature_wet_bulb": {"method": None},
            "mol_frac_gas_phase": {"method": None},
            "mole_frac_1_gas_phase": {"method": None},
            "mole_frac_1_liq_phase": {"method": None},
        })

        obj.add_default_units({
            "time": units.s,
            "length": units.m,
            "mass": units.kg,
            "amount": units.mol,
            "temperature": units.K,
        })

@declare_process_block_class("HAirStateBlock", block_class=StateBlock)
class HAirStateBlockData(StateBlockData):
    def build(self):
        super().build()
        self._build_state_vars()
        self._load_surrogate()
        self._build_property_expressions()

    def _build_state_vars(self):
        self.flow_mol = Var(domain=NonNegativeReals, initialize=1.0, units=units.mol / units.s)
        self.temperature = Var(domain=NonNegativeReals, initialize=300, units=units.K)
        self.pressure = Var(domain=NonNegativeReals, initialize=101325, units=units.Pa)
        self.mole_frac_comp = Var(self.params.component_list, domain=Reals, bounds=(0, 1), initialize=0.5, units=units.dimensionless)

    def _load_surrogate(self):
        self.relative_humidity = Var(initialize=0.5, units=units.dimensionless)
        self.temperature_wet_bulb = Var(initialize=290, units=units.K)

        self.enth_mol_phase = Var(self.params.phase_list, initialize=15000, units=units.J / units.mol)
        self.entr_mol_phase = Var(self.params.phase_list, initialize=40, units=units.J / units.mol / units.K)
        self.vol_mol_phase = Var(self.params.phase_list, initialize=1e-3, units=units.m**3 / units.mol)

        self.mol_frac_gas_phase = Var(initialize=0.8, units=units.dimensionless)
        self.mole_frac_1_gas_phase = Var(initialize=0.02, units=units.dimensionless)
        self.mole_frac_1_liq_phase = Var(initialize=0.98, units=units.dimensionless)

        script_dir = os.path.dirname(__file__)
        path = os.path.abspath(os.path.join(script_dir, "../results/pysmo_humid_air_props_surrogate.json"))
        if not os.path.exists(path):
            warnings.warn("Surrogate file not found.")
            return

        self.pysmo_surrogate = PysmoSurrogate.load_from_file(path)
        self.surrogate = SurrogateBlock()
        self.surrogate.build_model(
            self.pysmo_surrogate,
            input_vars=[
                self.temperature,
                self.pressure,
                self.mole_frac_comp['water']
            ],
            output_vars=[
                self.relative_humidity,
                self.temperature_wet_bulb,
                self.enth_mol_phase['Vap'],
                self.entr_mol_phase['Vap'],
                self.vol_mol_phase['Vap'],
                self.enth_mol_phase['Liq'],
                self.entr_mol_phase['Liq'],
                self.vol_mol_phase['Liq'],
                self.mol_frac_gas_phase,
                self.mole_frac_1_gas_phase,
                self.mole_frac_1_liq_phase
            ]
        )

    def _build_property_expressions(self):
        self.enth_mol = Expression(expr=(
            self.mol_frac_gas_phase * self.enth_mol_phase['Vap']
            + (1 - self.mol_frac_gas_phase) * self.enth_mol_phase['Liq']
        ))

        self.entr_mol = Expression(expr=(
            self.mol_frac_gas_phase * self.entr_mol_phase['Vap']
            + (1 - self.mol_frac_gas_phase) * self.entr_mol_phase['Liq']
        ))

        self.vol_mol = Expression(expr=(
            self.mol_frac_gas_phase * self.vol_mol_phase['Vap']
            + (1 - self.mol_frac_gas_phase) * self.vol_mol_phase['Liq']
        ))

    def define_state_vars(self):
        return {
            "flow_mol": self.flow_mol,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "mole_frac_comp": self.mole_frac_comp,
        }

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_terms(self, p, j):
        return self.flow_mol * self.mole_frac_comp[j]

    def get_enthalpy_flow_terms(self, p):
        return self.flow_mol * self.enth_mol
