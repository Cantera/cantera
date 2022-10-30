from contextlib import nullcontext
from dataclasses import dataclass
from typing import Optional, Tuple, Dict
import sys

import pytest
import cantera.with_units as ctu
import cantera as ct
try:
    from pint.testing import assert_allclose
except ModuleNotFoundError:
    # pint.testing was introduced in pint 0.20
    from pint.testsuite.helpers import assert_quantity_almost_equal as assert_allclose


@pytest.fixture(scope="function")
def ideal_gas():
    return ctu.Solution("h2o2.yaml")


@pytest.fixture(scope="function")
def pure_fluid():
    return ctu.Water()


@pytest.fixture(params=["ideal_gas", "pure_fluid"])
def generic_phase(request):
    return request.getfixturevalue(request.param)


def test_setting_basis_units_fails(generic_phase):
    # Python 3.10 includes the name of the attribute which was improperly used as a
    # setter. Earlier versions have just a generic error message.
    if sys.version_info.minor < 10:
        match = "set attribute"
    else:
        match = "basis_units"
    with pytest.raises(AttributeError, match=match):
        generic_phase.basis_units = "some random string"


def test_mass_basis(generic_phase):
    """Check that mass basis units have kg and the generic getter returns the same
    value as the mass-specific getter."""
    generic_phase.basis = "mass"
    assert generic_phase.basis_units == "kg"
    assert_allclose(generic_phase.density_mass, generic_phase.density)
    assert_allclose(generic_phase.enthalpy_mass, generic_phase.h)
    assert_allclose(generic_phase.entropy_mass, generic_phase.s)
    assert_allclose(generic_phase.int_energy_mass, generic_phase.u)
    assert_allclose(generic_phase.volume_mass, generic_phase.v)
    assert_allclose(generic_phase.gibbs_mass, generic_phase.g)
    assert_allclose(generic_phase.cp_mass, generic_phase.cp)
    assert_allclose(generic_phase.cv_mass, generic_phase.cv)


def test_molar_basis(generic_phase):
    """Check that molar basis units have kmol and the generic getter returns the
    same value as the molar-specific getter."""
    generic_phase.basis = "molar"
    assert generic_phase.basis_units == "kmol"
    assert_allclose(generic_phase.density_mole, generic_phase.density)
    assert_allclose(generic_phase.enthalpy_mole, generic_phase.h)
    assert_allclose(generic_phase.entropy_mole, generic_phase.s)
    assert_allclose(generic_phase.int_energy_mole, generic_phase.u)
    assert_allclose(generic_phase.volume_mole, generic_phase.v)
    assert_allclose(generic_phase.gibbs_mole, generic_phase.g)
    assert_allclose(generic_phase.cp_mole, generic_phase.cp)
    assert_allclose(generic_phase.cv_mole, generic_phase.cv)


@dataclass(frozen=True)
class Dimensions:
    name: str
    mass: Optional[float] = None
    length: Optional[float] = None
    time: Optional[float] = None
    substance: Optional[float] = None
    temperature: Optional[float] = None
    current: Optional[float] = None
    dimensions: Tuple[str, ...] = (
        "mass",
        "length",
        "time",
        "substance",
        "temperature",
        "current",
    )

    def __str__(self):
        return self.name

    def inverse(self, name: Optional[str] = None) -> "Dimensions":
        dimensionality = {}
        for dimension in self.dimensions:
            value = getattr(self, dimension)
            if value is not None:
                dimensionality[dimension] = -1 * value
        new_name = name if name is not None else f"inverse_{self.name}"
        return Dimensions(name=new_name, **dimensionality)

    def as_dict(self) -> Dict[str, float]:
        """Add the square brackets around the dimension for comparison with pint"""
        dimensionality = {}
        for dimension in self.dimensions:
            value = getattr(self, dimension)
            if value is not None:
                dimensionality[f"[{dimension}]"] = value
        return dimensionality


temperature = Dimensions(temperature=1, name="temperature")
pressure = Dimensions(mass=1, length=-1, time=-2, name="pressure")
isothermal_compressiblity = pressure.inverse()
inverse_temperature = temperature.inverse()
atomic_molecular_weights = Dimensions("atomic_molecular_weights", mass=1, substance=-1)
mole_mass_fractions = Dimensions(name="mole_mass_fractions")
chemical_potential = Dimensions(
    name="chemical_potential", mass=1, length=2, time=-2, substance=-1
)
electric_potential = Dimensions(
    name="electric_potential", mass=1, length=2, time=-3, current=-1
)
concentrations_like = Dimensions(name="concentrations_like", substance=1, length=-3)
molar_volume = Dimensions(name="volume_mole", substance=-1, length=3)
volume_mass = Dimensions(name="volume_mass", mass=-1, length=3)
density_mass = volume_mass.inverse(name="density_mass")
mass_basis_energy_like = Dimensions(name="mass_basis_energy_like", length=2, time=-2)
mass_basis_entropy_like = Dimensions(
    name="mass_basis_entropy_like", length=2, time=-2, temperature=-1
)
molar_basis_energy_like = Dimensions(
    name="molar_basis_energy_like", mass=1, length=2, time=-2, substance=-1
)
molar_basis_entropy_like = Dimensions(
    name="molar_basis_entropy_like",
    mass=1,
    length=2,
    time=-2,
    substance=-1,
    temperature=-1,
)


def yield_dimensions():
    """Yield pytest.param instances with the dimensions"""
    # This dictionary maps the dimensions to the relevant property names
    dims: Dict[Dimensions, Tuple[str, ...]] = {
        temperature: ("T", "max_temp", "min_temp"),
        inverse_temperature: ("thermal_expansion_coeff",),
        pressure: ("P", "reference_pressure"),
        isothermal_compressiblity: ("isothermal_compressibility",),
        atomic_molecular_weights: (
            "atomic_weight",
            "molecular_weights",
            "mean_molecular_weight",
        ),
        mole_mass_fractions: ("X", "Y"),
        chemical_potential: ("chemical_potentials", "electrochemical_potentials"),
        electric_potential: ("electric_potential",),
        concentrations_like: ("concentrations", "density_mole"),
        molar_volume: ("volume_mole", "partial_molar_volumes"),
        volume_mass: ("volume_mass",),
        density_mass: ("density_mass",),
        mass_basis_energy_like: ("enthalpy_mass", "int_energy_mass", "gibbs_mass"),
        mass_basis_entropy_like: ("entropy_mass", "cp_mass", "cv_mass"),
        molar_basis_energy_like: (
            "enthalpy_mole",
            "int_energy_mole",
            "gibbs_mole",
            "partial_molar_enthalpies",
            "partial_molar_int_energies",
        ),
        molar_basis_entropy_like: (
            "entropy_mole",
            "cp_mole",
            "cv_mole",
            "partial_molar_cp",
            "partial_molar_entropies",
        ),
    }
    for dimension, props in dims.items():
        for prop in props:
            yield pytest.param(prop, dimension.as_dict(), id=f"{dimension}-{prop}")


@pytest.mark.parametrize("prop,dimensions", yield_dimensions())
def test_dimensions(generic_phase, prop, dimensions):
    pint_dim = dict(getattr(generic_phase, prop).dimensionality)
    assert pint_dim == dimensions


@pytest.mark.parametrize(
    "phase",
    (
        pytest.param(ctu.Solution("liquidvapor.yaml", "heptane"), id="Solution"),
        pytest.param(ctu.Water(), id="PureFluid"),
    ),
)
def test_purefluid_dimensions(phase):
    # Test some dimensions that weren't tested as part of the Solution tests
    # Create and test a liquidvapor phase in a Solution object, since an ideal gas phase
    # doesn't implement saturation or critical properties.
    assert dict(phase.T_sat.dimensionality) == temperature.as_dict()
    assert dict(phase.critical_temperature.dimensionality) == temperature.as_dict()
    assert dict(phase.P_sat.dimensionality) == pressure.as_dict()
    assert dict(phase.critical_pressure.dimensionality) == pressure.as_dict()
    assert dict(phase.critical_density.dimensionality) == density_mass.as_dict()


def yield_prop_pairs():
    pairs = [
        pytest.param(("TP", "T", "P"), id="TP"),
        pytest.param(("SP", "s", "P"), id="SP"),
        pytest.param(("UV", "u", "v"), id="UV"),
        pytest.param(("DP", "density", "P"), id="DP"),
        pytest.param(("HP", "h", "P"), id="HP"),
        pytest.param(("SV", "s", "v"), id="SV"),
        pytest.param(("TD", "T", "density"), id="TD"),
    ]
    yield from pairs


def yield_prop_triples():
    for pair in yield_prop_pairs():
        values = pair.values[0]
        yield pytest.param(
            (values[0] + "X", *values[1:], "X"),
            id=pair.id + "X",
        )
        yield pytest.param(
            (values[0] + "Y", *values[1:], "Y"),
            id=pair.id + "Y",
        )


def yield_prop_pairs_and_triples():
    yield from yield_prop_pairs()
    yield from yield_prop_triples()


@pytest.fixture
def TD_in_the_right_basis(request):
    if request.param == "mass":
        T = ctu.Q_(500, "K")
        rho = ctu.Q_(1.5, "kg/m**3")
    elif request.param == "molar":
        T = ctu.Q_(750, "K")
        rho = ctu.Q_(0.02, "kmol/m**3")
    return (T, rho)


@pytest.fixture
def initial_TDY(request, generic_phase):
    generic_phase.basis = request.param
    return generic_phase.TDY


@pytest.fixture
def some_setters_arent_implemented_for_purefluid(request):
    pair_or_triple = request.getfixturevalue("props")[0]
    is_pure_fluid = isinstance(request.getfixturevalue("generic_phase"), ctu.PureFluid)
    if is_pure_fluid and pair_or_triple.startswith("DP"):
        request.applymarker(
            pytest.mark.xfail(
                raises=ctu.CanteraError,
                reason=f"The {pair_or_triple} method isn't implemented",
            )
        )


# The parameterization is done here with the indirect kwarg to make sure that the same
# basis is passed to both fixtures. The alternative is to use the params kwarg to the
# fixture decorator, which would give us (mass, molar) basis pairs, and that doesn't
# make sense.
@pytest.mark.parametrize(
    "TD_in_the_right_basis,initial_TDY",
    [
        pytest.param("mass", "mass", id="mass"),
        pytest.param("molar", "molar", id="molar"),
    ],
    indirect=True,
)
@pytest.mark.parametrize("props", yield_prop_pairs_and_triples())
@pytest.mark.usefixtures("some_setters_arent_implemented_for_purefluid")
def test_setters(generic_phase, TD_in_the_right_basis, initial_TDY, props):
    pair_or_triple = props[0]
    if isinstance(generic_phase, ctu.PureFluid):
        Y_1 = ctu.Q_([1.0], "dimensionless")
    else:
        Y_1 = ctu.Q_(
            [0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2, 0.0], "dimensionless"
        )
    generic_phase.TDY = *TD_in_the_right_basis, Y_1

    # Use TDY setting to get the properties at the modified state
    new_props = getattr(generic_phase, pair_or_triple)

    # Reset to the initial state so that the next state setting actually has to do
    # something.
    generic_phase.TDY = initial_TDY

    # If we're only setting a pair of properties, reset the mass fractions to the
    # expected state before using the pair to set.
    if len(pair_or_triple) == 2:
        generic_phase.Y = Y_1

    # Use the test pair or triple to set the state and assert that the
    # natural properties are equal to the modified state
    setattr(generic_phase, pair_or_triple, new_props)
    T_1, rho_1 = TD_in_the_right_basis
    assert_allclose(generic_phase.T, T_1)
    assert_allclose(generic_phase.density, rho_1)
    assert_allclose(generic_phase.Y, Y_1)


@pytest.mark.parametrize("props", yield_prop_triples())
@pytest.mark.usefixtures("some_setters_arent_implemented_for_purefluid")
def test_setters_hold_constant(generic_phase, props):
    triple, first, second, third = props

    # Set an arbitrary initial state
    if generic_phase.n_species == 1:
        generic_phase.X = ctu.Q_([1.0], "dimensionless")
        composition = ctu.Q_([1.0], "dimensionless")
    else:
        generic_phase.X = "H2O:0.1, O2:0.95, AR:3.0"
        composition = "H2:0.1, O2:1.0, AR:3.0"
    generic_phase.TD = ctu.Q_(1000, "K"), ctu.Q_(1.5, "kg/m**3")
    property_3 = getattr(generic_phase, third)

    # Change to another arbitrary state and store values to compare when a property
    # isn't changed
    reset_state = (ctu.Q_(500, "K"), ctu.Q_(2.5, "kg/m**3"), composition)
    generic_phase.TDX = reset_state
    first_val, second_val, third_val = getattr(generic_phase, triple)

    setattr(generic_phase, triple, (None, None, property_3))
    assert_allclose(getattr(generic_phase, first), first_val)
    assert_allclose(getattr(generic_phase, second), second_val)
    assert_allclose(getattr(generic_phase, third), property_3)

    generic_phase.TDX = reset_state
    setattr(generic_phase, triple, (None, None, None))
    assert_allclose(getattr(generic_phase, first), first_val)
    assert_allclose(getattr(generic_phase, second), second_val)
    assert_allclose(getattr(generic_phase, third), third_val)


@pytest.mark.parametrize("props", yield_prop_pairs_and_triples())
@pytest.mark.parametrize(
    "basis,rho_0",
    [
        pytest.param("mass", ctu.Q_(0.7, "kg/m**3"), id="mass"),
        pytest.param("molar", ctu.Q_(0.01, "kmol/m**3"), id="molar"),
    ],
)
def test_multi_prop_getters_are_equal_to_single(generic_phase, props, basis, rho_0):
    pair_or_triple, first, second, *third = props
    generic_phase.basis = basis
    if generic_phase.n_species != 1:
        generic_phase.Y = "H2:0.1, H2O2:0.1, AR:0.8"
    generic_phase.TD = ctu.Q_(350.0, "K"), rho_0
    first_value, second_value, *third_value = getattr(generic_phase, pair_or_triple)
    assert_allclose(getattr(generic_phase, first), first_value)
    assert_allclose(getattr(generic_phase, second), second_value)
    if third:
        assert_allclose(getattr(generic_phase, third[0]), third_value[0])


@pytest.mark.parametrize("pair", yield_prop_pairs())
def test_set_pair_without_units_is_an_error(generic_phase, pair):
    with pytest.raises(ctu.CanteraError, match="an instance of a pint"):
        setattr(generic_phase, pair[0], [300, None])


@pytest.fixture
def third_prop_sometimes_fails(request):
    # Setting these property triples with None, None, [all equal fractions] results in
    # a convergence error in setting the state. We don't care that it fails because
    # this is just testing what happens to the value that gets passed in to the
    # with_units setters in terms of conversion to/from values that don't have units
    # already attached to them.
    prop_triple = request.getfixturevalue("triple")[0]
    is_purefluid = isinstance(request.getfixturevalue("generic_phase"), ctu.PureFluid)
    if prop_triple in ("SPX", "UVX", "UVY", "HPX", "HPY", "SVX") and not is_purefluid:
        return pytest.raises(ctu.CanteraError, match="ThermoPhase::setState")
    elif prop_triple in ("DPX", "DPY") and is_purefluid:
        return pytest.raises(ctu.CanteraError, match="setState_RP")
    else:
        return nullcontext()


@pytest.mark.parametrize("triple", yield_prop_triples())
def test_set_triple_without_units_is_an_error(
    generic_phase,
    triple,
    third_prop_sometimes_fails,
):
    value_1 = [300, None, [1] * generic_phase.n_species]
    with pytest.raises(ctu.CanteraError, match="an instance of a pint"):
        setattr(generic_phase, triple[0], value_1)

    value_3 = [None, None, [1] * generic_phase.n_species]
    with third_prop_sometimes_fails:
        setattr(generic_phase, triple[0], value_3)


@pytest.fixture
def pure_fluid_in_vapordome(pure_fluid):
    pure_fluid.TQ = None, ctu.Q_(0.5, "dimensionless")
    return pure_fluid


def test_set_Q(pure_fluid_in_vapordome):
    P = pure_fluid_in_vapordome.P
    T = pure_fluid_in_vapordome.T
    pure_fluid_in_vapordome.Q = ctu.Q_(0.6, "dimensionless")
    assert_allclose(pure_fluid_in_vapordome.Q, 0.6 * ctu.units.dimensionless)
    assert_allclose(pure_fluid_in_vapordome.T, T)
    assert_allclose(pure_fluid_in_vapordome.P, P)

    pure_fluid_in_vapordome.Q = None
    assert_allclose(pure_fluid_in_vapordome.Q, 0.6 * ctu.units.dimensionless)
    assert_allclose(pure_fluid_in_vapordome.T, T)
    assert_allclose(pure_fluid_in_vapordome.P, P)

    with pytest.raises(ctu.CanteraError, match="an instance of a pint"):
        pure_fluid_in_vapordome.Q = 0.5


def yield_purefluid_only_setters():
    props = [
        pytest.param(("TPQ", "T", "P", "Q"), id="TPQ"),
        pytest.param(("TQ", "T", "Q"), id="TQ"),
        pytest.param(("PQ", "P", "Q"), id="PQ"),
        pytest.param(("PV", "P", "v"), id="PV"),
        pytest.param(("SH", "s", "h"), id="SH"),
        pytest.param(("ST", "s", "T"), id="ST"),
        pytest.param(("TH", "T", "h"), id="TH"),
        pytest.param(("TV", "T", "v"), id="TV"),
        pytest.param(("VH", "v", "h"), id="VH"),
        pytest.param(("UP", "u", "P"), id="UP"),
    ]
    yield from props


def yield_purefluid_only_getters():
    props = [
        pytest.param(("DPQ", "density", "P", "Q"), id="DPQ"),
        pytest.param(("HPQ", "h", "P", "Q"), id="HPQ"),
        pytest.param(("SPQ", "s", "P", "Q"), id="SPQ"),
        pytest.param(("SVQ", "s", "v", "Q"), id="SVQ"),
        pytest.param(("TDQ", "T", "density", "Q"), id="TDQ"),
        pytest.param(("UVQ", "u", "v", "Q"), id="UVQ"),
    ]
    yield from props


def yield_all_purefluid_only_props():
    yield from yield_purefluid_only_setters()
    yield from yield_purefluid_only_getters()


@pytest.mark.parametrize("prop", yield_purefluid_only_setters())
def test_set_without_units_is_error_purefluid(prop, pure_fluid):
    value = [None] * (len(prop[0]) - 1) + [0.5]
    with pytest.raises(ctu.CanteraError, match="an instance of a pint"):
        # Don't use append here, because append returns None which would be
        # passed to setattr
        setattr(pure_fluid, prop[0], value)


@pytest.mark.parametrize("props", yield_all_purefluid_only_props())
def test_multi_prop_getters_purefluid(pure_fluid, props):
    pair_or_triple, first, second, *third = props
    first_value, second_value, *third_value = getattr(pure_fluid, pair_or_triple)
    assert_allclose(getattr(pure_fluid, first), first_value)
    assert_allclose(getattr(pure_fluid, second), second_value)
    if third:
        assert_allclose(getattr(pure_fluid, third[0]), third_value[0])


@pytest.mark.parametrize("props", yield_purefluid_only_setters())
def test_setters_purefluid(props, pure_fluid):
    # Only need to run this for a single pure fluid
    initial_TD = pure_fluid.TD
    pair_or_triple = props[0]

    T_1 = ctu.Q_(500, "K")
    if pair_or_triple in ("SH", "TH"):
        # This state is able to converge for these setters, whereas the state below
        # does not converge
        rho_1 = ctu.Q_(1000, "kg/m**3")
    else:
        # This state is located inside the vapor dome to be able to test the
        # TQ and PQ setters
        rho_1 = ctu.Q_(25.93245092697775, "kg/m**3")

    # Use TD setting to get the properties at the modified state
    pure_fluid.TD = T_1, rho_1
    new_props = getattr(pure_fluid, pair_or_triple)
    print(new_props)

    # Reset to the initial state so that the next state setting actually has to do
    # something.
    pure_fluid.TD = initial_TD

    # Use the test pair or triple to set the state and assert that the
    # natural properties are equal to the modified state
    setattr(pure_fluid, pair_or_triple, new_props)
    assert_allclose(pure_fluid.T, T_1)
    assert_allclose(pure_fluid.density, rho_1)


@pytest.mark.parametrize("prop", ("X", "Y"))
def test_X_Y_setters_with_none(generic_phase, prop):
    comparison = getattr(generic_phase, prop)
    setattr(generic_phase, prop, None)
    # Assert that the value hasn't changed
    assert_allclose(comparison, getattr(generic_phase, prop))


@pytest.mark.parametrize("prop", ("X", "Y"))
def test_X_Y_setters_without_units_works(generic_phase, prop):
    composition = f"{generic_phase.species_names[0]}:1"
    setattr(generic_phase, prop, composition)
    assert_allclose(getattr(generic_phase, prop)[0], ctu.Q_([1], "dimensionless"))


def test_thermophase_properties_exist(ideal_gas):
    # Since the Solution class in the with_units subpackage only implements
    # the ThermoPhase interface for now, instantiate a regular ThermoPhase
    # to compare the attributes and make sure all of them exist on the with_units
    # object
    tp = ct.ThermoPhase("h2o2.yaml")

    for attr in dir(tp):
        if attr.startswith("_"):  # or attr in skip:
            continue

        try:
            getattr(tp, attr)
        except (ct.CanteraError, ct.ThermoModelMethodError):
            continue

        assert hasattr(ideal_gas, attr)


def test_purefluid_properties_exist(pure_fluid):
    pf = ct.PureFluid("liquidvapor.yaml", "water")
    for attr in dir(pf):
        if attr.startswith("_"):
            continue

        try:
            getattr(pf, attr)
        except (ct.CanteraError, ct.ThermoModelMethodError):
            continue

        assert hasattr(pure_fluid, attr)
