# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

import sys as _sys
import importlib.abc as _abc
import importlib.machinery as _machinery
import importlib.util as _util

# Chooses the right init function
# See https://stackoverflow.com/a/52714500
class _CythonPackageMetaPathFinder(_abc.MetaPathFinder):
    def __init__(self, name_filter):
        super().__init__()
        self.name_filter = name_filter

    def find_spec(self, fullname, path, target=None):
        if fullname.startswith(self.name_filter):
            # use this extension-file but PyInit-function of another module:
            loader = _machinery.ExtensionFileLoader(fullname, __file__)
            return _util.spec_from_loader(fullname, loader)


# Inject custom finder/loaders into sys.meta_path:
_sys.meta_path.append(_CythonPackageMetaPathFinder("cantera."))

# Import the contents of the individual .pyx files
from ._onedim import (
    AxisymmetricFlow, Boundary1D, Domain1D, FlowBase, FreeFlow, Inlet1D, Outlet1D,
    OutletReservoir1D, ReactingSurface1D, Sim1D, Surface1D, SymmetryPlane1D,
    UnstrainedFlow,
)
from ._utils import (
    AnyMap, CanteraError,
    __git_commit__, __sundials_version__, __version__,
    add_directory, appdelete, debug_mode_enabled, get_data_directories, hdf_support,
    make_deprecation_warnings_fatal, print_stack_trace_on_segfault,
    suppress_deprecation_warnings, suppress_thermo_warnings, use_legacy_rate_constants,
    use_sparse,
)
from .constants import (
    avogadro, boltzmann, electron_charge, electron_mass, epsilon_0, faraday,
    gas_constant, light_speed, one_atm, permeability_0, planck, stefan_boltzmann,
)
from .delegator import extension
from .func1 import Func1, Tabulated1
from .jacobians import (
    AdaptivePreconditioner, BandedJacobian, EigenSparseDirectJacobian,
    EigenSparseJacobian, SystemJacobian,
)
from .kinetics import Kinetics, InterfaceKinetics
from .mixture import Mixture
from .reaction import (
    Arrhenius, ArrheniusRate, ArrheniusRateBase, BlowersMaselRate, ChebyshevRate,
    CustomRate, ElectronCollisionPlasmaRate, ExtensibleRate, ExtensibleRateData,
    FalloffRate, InterfaceArrheniusRate, InterfaceBlowersMaselRate,
    InterfaceRateBase, LinearBurkeRate, LindemannRate, PlogRate, Reaction,
    ReactionRate, SriRate, StickRateBase, StickingArrheniusRate,
    StickingBlowersMaselRate, ThirdBody, TroeRate, TsangRate, TwoTempPlasmaRate,
)
from .reactionpath import ReactionPathDiagram
from .reactor import (
    ConnectorNode, ConstPressureMoleReactor, ConstPressureReactor,
    ExtensibleConstPressureMoleReactor, ExtensibleConstPressureReactor,
    ExtensibleIdealGasConstPressureMoleReactor,
    ExtensibleIdealGasConstPressureReactor, ExtensibleIdealGasMoleReactor,
    ExtensibleIdealGasReactor, ExtensibleMoleReactor, ExtensibleReactor,
    FlowDevice, FlowReactor, IdealGasConstPressureMoleReactor,
    IdealGasConstPressureReactor, IdealGasMoleReactor, IdealGasReactor,
    MassFlowController, MoleReactor, PressureController, Reactor, ReactorBase,
    ReactorNet, ReactorSurface, Reservoir, Valve, Wall, WallBase,
)
from .solutionbase import SolutionArrayBase
from .speciesthermo import (
    ConstantCp, Mu0Poly, Nasa9PolyMultiTempRegion, NasaPoly2, ShomatePoly2,
    SpeciesThermo,
)
from .thermo import (
    Element, InterfacePhase, PureFluid, Species, ThermoModelMethodError, ThermoPhase,
)
from .transport import DustyGasTransport, GasTransportData, Transport
from .units import UnitStack, UnitSystem, Units
from .yamlwriter import YamlWriter


# Custom finder/loader no longer needed, so remove it
_sys.meta_path.pop()
