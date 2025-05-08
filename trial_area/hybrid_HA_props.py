# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Param,
    Reals,
    Set,
    value,
    Var,
    NonNegativeReals,
    units,
)
from pyomo.opt import SolverFactory, TerminationCondition

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    LiquidPhase,
    Component,
)
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.misc import extract_data
from idaes.core.solvers import get_solver
from pyomo.util.check_units import assert_units_consistent
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate

from pyomo.util.model_size import build_model_size_report

# Some more information about this module
__author__ = Ben

# Set up logger
_log = logging.getLogger(__name__)