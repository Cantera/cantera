# Source files for the Cantera library
# Organized by subdirectory

# Base module
base_sources = files(
  'base/AnyMap.cpp',
  'base/Array.cpp',
  'base/ExtensionManager.cpp',
  'base/ExtensionManagerFactory.cpp',
  'base/Interface.cpp',
  'base/Solution.cpp',
  'base/SolutionArray.cpp',
  'base/Storage.cpp',
  'base/Units.cpp',
  'base/ValueCache.cpp',
  'base/YamlWriter.cpp',
  'base/application.cpp',
  'base/checkFinite.cpp',
  'base/clockWC.cpp',
  'base/ctexceptions.cpp',
  'base/global.cpp',
  'base/stringUtils.cpp',
)

# Thermo module
thermo_sources = files(
  'thermo/BinarySolutionTabulatedThermo.cpp',
  'thermo/ConstCpPoly.cpp',
  'thermo/CoverageDependentSurfPhase.cpp',
  'thermo/DebyeHuckel.cpp',
  'thermo/EEDFTwoTermApproximation.cpp',
  'thermo/Elements.cpp',
  'thermo/GibbsExcessVPSSTP.cpp',
  'thermo/HMWSoln.cpp',
  'thermo/IdealGasPhase.cpp',
  'thermo/IdealMolalSoln.cpp',
  'thermo/IdealSolidSolnPhase.cpp',
  'thermo/IdealSolnGasVPSS.cpp',
  'thermo/MargulesVPSSTP.cpp',
  'thermo/MixtureFugacityTP.cpp',
  'thermo/MolalityVPSSTP.cpp',
  'thermo/Mu0Poly.cpp',
  'thermo/MultiSpeciesThermo.cpp',
  'thermo/Nasa9Poly1.cpp',
  'thermo/Nasa9PolyMultiTempRegion.cpp',
  'thermo/NasaPoly2.cpp',
  'thermo/PDSS.cpp',
  'thermo/PDSSFactory.cpp',
  'thermo/PDSS_ConstVol.cpp',
  'thermo/PDSS_HKFT.cpp',
  'thermo/PDSS_SSVol.cpp',
  'thermo/PDSS_Water.cpp',
  'thermo/PengRobinson.cpp',
  'thermo/Phase.cpp',
  'thermo/PlasmaPhase.cpp',
  'thermo/PureFluidPhase.cpp',
  'thermo/RedlichKisterVPSSTP.cpp',
  'thermo/RedlichKwongMFTP.cpp',
  'thermo/SingleSpeciesTP.cpp',
  'thermo/Species.cpp',
  'thermo/SpeciesThermoFactory.cpp',
  'thermo/SpeciesThermoInterpType.cpp',
  'thermo/StoichSubstance.cpp',
  'thermo/SurfPhase.cpp',
  'thermo/ThermoFactory.cpp',
  'thermo/ThermoPhase.cpp',
  'thermo/VPStandardStateTP.cpp',
  'thermo/WaterProps.cpp',
  'thermo/WaterPropsIAPWS.cpp',
  'thermo/WaterPropsIAPWSphi.cpp',
  'thermo/WaterSSTP.cpp',
)

# TPX module
tpx_sources = files(
  'tpx/CarbonDioxide.cpp',
  'tpx/HFC134a.cpp',
  'tpx/Heptane.cpp',
  'tpx/Hydrogen.cpp',
  'tpx/Methane.cpp',
  'tpx/Nitrogen.cpp',
  'tpx/Oxygen.cpp',
  'tpx/Sub.cpp',
  'tpx/Water.cpp',
  'tpx/utils.cpp',
)

# Equilibrium module
equil_sources = files(
  'equil/BasisOptimize.cpp',
  'equil/ChemEquil.cpp',
  'equil/MultiPhase.cpp',
  'equil/MultiPhaseEquil.cpp',
  'equil/vcs_MultiPhaseEquil.cpp',
  'equil/vcs_VolPhase.cpp',
  'equil/vcs_solve.cpp',
  'equil/vcs_solve_TP.cpp',
  'equil/vcs_util.cpp',
)

# Numerics module
numerics_sources = files(
  'numerics/AdaptivePreconditioner.cpp',
  'numerics/BandMatrix.cpp',
  'numerics/CVodesIntegrator.cpp',
  'numerics/DenseMatrix.cpp',
  'numerics/EigenSparseDirectJacobian.cpp',
  'numerics/EigenSparseJacobian.cpp',
  'numerics/Func1.cpp',
  'numerics/Func1Factory.cpp',
  'numerics/FuncEval.cpp',
  'numerics/IdasIntegrator.cpp',
  'numerics/Integrators.cpp',
  'numerics/SteadyStateSystem.cpp',
  'numerics/SystemJacobianFactory.cpp',
  'numerics/funcs.cpp',
  'numerics/polyfit.cpp',
)

# Kinetics module
kinetics_sources = files(
  'kinetics/Arrhenius.cpp',
  'kinetics/BlowersMaselRate.cpp',
  'kinetics/BulkKinetics.cpp',
  'kinetics/ChebyshevRate.cpp',
  'kinetics/Custom.cpp',
  'kinetics/ElectronCollisionPlasmaRate.cpp',
  'kinetics/Falloff.cpp',
  'kinetics/Group.cpp',
  'kinetics/InterfaceKinetics.cpp',
  'kinetics/InterfaceRate.cpp',
  'kinetics/Kinetics.cpp',
  'kinetics/KineticsFactory.cpp',
  'kinetics/LinearBurkeRate.cpp',
  'kinetics/PlogRate.cpp',
  'kinetics/Reaction.cpp',
  'kinetics/ReactionPath.cpp',
  'kinetics/ReactionRateDelegator.cpp',
  'kinetics/ReactionRateFactory.cpp',
  'kinetics/TwoTempPlasmaRate.cpp',
)

# Transport module
transport_sources = files(
  'transport/DustyGasTransport.cpp',
  'transport/GasTransport.cpp',
  'transport/HighPressureGasTransport.cpp',
  'transport/IonGasTransport.cpp',
  'transport/MMCollisionInt.cpp',
  'transport/MixTransport.cpp',
  'transport/MultiTransport.cpp',
  'transport/Transport.cpp',
  'transport/TransportData.cpp',
  'transport/TransportFactory.cpp',
  'transport/WaterTransport.cpp',
)

# OneD module
oned_sources = files(
  'oneD/Boundary1D.cpp',
  'oneD/Domain1D.cpp',
  'oneD/DomainFactory.cpp',
  'oneD/Flow1D.cpp',
  'oneD/IonFlow.cpp',
  'oneD/MultiJac.cpp',
  'oneD/MultiNewton.cpp',
  'oneD/OneDim.cpp',
  'oneD/Sim1D.cpp',
  'oneD/refine.cpp',
)

# ZeroD module
zerod_sources = files(
  'zeroD/ConnectorFactory.cpp',
  'zeroD/ConnectorNode.cpp',
  'zeroD/ConstPressureMoleReactor.cpp',
  'zeroD/ConstPressureReactor.cpp',
  'zeroD/FlowDevice.cpp',
  'zeroD/FlowReactor.cpp',
  'zeroD/IdealGasConstPressureMoleReactor.cpp',
  'zeroD/IdealGasConstPressureReactor.cpp',
  'zeroD/IdealGasMoleReactor.cpp',
  'zeroD/IdealGasReactor.cpp',
  'zeroD/MoleReactor.cpp',
  'zeroD/Reactor.cpp',
  'zeroD/ReactorBase.cpp',
  'zeroD/ReactorFactory.cpp',
  'zeroD/ReactorNet.cpp',
  'zeroD/ReactorSurface.cpp',
  'zeroD/Wall.cpp',
  'zeroD/flowControllers.cpp',
)

# Extensions module (minimal for now)
extensions_sources = files(
  'extensions/PythonExtensionManager.cpp',
  'extensions/canteraShared.cpp',
  'extensions/pythonShim.cpp',
)

# Combine all sources
cantera_sources = [
  base_sources,
  thermo_sources,
  tpx_sources,
  equil_sources,
  numerics_sources,
  kinetics_sources,
  transport_sources,
  oned_sources,
  zerod_sources,
  extensions_sources,
]

# Additional compile definitions for specific files
application_cpp_args = ['-DCANTERA_DATA="' + datadir + '"']
global_cpp_args = ['-DGIT_COMMIT="' + git_commit + '"']

# Build application.cpp with special define
application_lib = static_library('cantera_application',
  'base/application.cpp',
  include_directories: cantera_include_dirs,
  cpp_args: application_cpp_args,
  dependencies: [
    eigen_dep,
    fmt_dep,
    yamlcpp_dep,
    boost_dep,
    thread_dep,
    m_dep,
  ],
  pic: true,
)

# Build global.cpp with special define
global_lib = static_library('cantera_global',
  'base/global.cpp',
  include_directories: cantera_include_dirs,
  cpp_args: global_cpp_args,
  dependencies: [
    eigen_dep,
    fmt_dep,
    yamlcpp_dep,
    boost_dep,
    thread_dep,
    m_dep,
  ],
  pic: true,
)

# Collect dependencies
cantera_deps = [
  eigen_dep,
  fmt_dep,
  yamlcpp_dep,
  boost_dep,
  thread_dep,
  m_dep,
]

if sundials_dep.found()
  cantera_deps += [sundials_dep, sundials_idas_dep, sundials_nvecserial_dep]
endif

if use_lapack
  cantera_deps += [lapack_dep]
  if blas_dep.found()
    cantera_deps += [blas_dep]
  endif
endif

if use_hdf5
  cantera_deps += [hdf5_dep, highfive_dep]
endif

# Build the shared library
cantera_lib = shared_library('cantera',
  cantera_sources,
  include_directories: cantera_include_dirs,
  dependencies: cantera_deps,
  link_whole: [application_lib, global_lib],
  install: true,
  version: meson.project_version(),
  soversion: version_major,
)

# Install headers
install_subdir('include/cantera',
  install_dir: get_option('includedir'),
  strip_directory: false,
)

# Declare dependency for use by other subprojects
cantera_dep = declare_dependency(
  link_with: cantera_lib,
  include_directories: cantera_include_dirs,
  dependencies: cantera_deps,
)

# pkg-config file
pkg = import('pkgconfig')
pkg.generate(cantera_lib,
  name: 'Cantera',
  description: 'Chemical kinetics, thermodynamics, and transport library',
  url: 'https://cantera.org',
  subdirs: 'cantera',
)
