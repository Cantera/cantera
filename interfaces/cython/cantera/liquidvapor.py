# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from . import PureFluid, _cantera


def Water():
    """
    Create a `PureFluid` object using the equation of state for water and the
    `WaterTransport` class for viscosity and thermal conductivity.

    The object returned by this method implements an accurate equation of
    state for water that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram. The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances.* Stanford: Stanford
    University, 1979. Print.

    whereas formulas for transport are taken from

    J. V. Sengers, J. T. R. Watson, *Improved International Formulations for
    the Viscosity and Thermal Conductivity of Water Substance,* J. Phys. Chem.
    Ref. Data, 15, 1291, 1986.

    For more details, see classes Cantera::PureFluid, tpx::water and
    Cantera::WaterTransport in the Cantera C++ source code documentation.
    """
    class WaterWithTransport(PureFluid, _cantera.Transport):
        __slots__ = ()

    return WaterWithTransport('liquidvapor.yaml', 'water',
                              transport_model='Water')


def LiquidWater():
    """
    Create a `PureFluid` object using the IAPWS Formulation 1995 for
    thermodynamic properties and the `WaterTransport` class for viscosity
    and thermal conductivity.

    The object returned by this method implements an accurate equation of
    state for water that can be used in the liquid and supercritical regions
    of the phase diagram. The equation of state is taken from

    W. Wagner, A. Pruss, *The IAPWS Formulation 1995 for the Thermodynamic
    Properties of Ordinary Water Substance for General and Scientific Use,*
    J. Phys. Chem. Ref. Dat, 31, 387, 2002.

    whereas formulas for transport are taken from

    J. V. Sengers, J. T. R. Watson, *Improved International Formulations for
    the Viscosity and Thermal Conductivity of Water Substance,* J. Phys. Chem.
    Ref. Data, 15, 1291, 1986.

    For more details, see classes Cantera::PureFluid, Cantera::WaterSSTP and
    Cantera::WaterTransport in the Cantera C++ source code documentation.
    """
    class WaterWithTransport(PureFluid, _cantera.Transport):
        __slots__ = ()

    return WaterWithTransport('liquidvapor.yaml', 'liquid-water',
                              transport_model='Water')


def Nitrogen():
    """
    Create a `PureFluid` object using the equation of state for nitrogen.

    The object returned by this method implements an accurate equation of
    state for nitrogen that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram. The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::nitrogen in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'nitrogen')


def Methane():
    """
    Create a `PureFluid` object using the equation of state for methane.

    The object returned by this method implements an accurate equation of
    state for methane that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram.  The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::methane in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'methane')


def Hydrogen():
    """
    Create a `PureFluid` object using the equation of state for hydrogen.

    The object returned by this method implements an accurate equation of
    state for hydrogen that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram. The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::hydrogen in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'hydrogen')


def Oxygen():
    """
    Create a `PureFluid` object using the equation of state for oxygen.

    The object returned by this method implements an accurate equation of
    state for oxygen that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram.  The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::oxygen in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'oxygen')


def Hfc134a():
    """
    Create a `PureFluid` object using the equation of state for HFC-134a.

    The object returned by this method implements an accurate equation of
    state for refrigerant HFC134a (R134a) that can be used in the liquid,
    vapor, saturated liquid/vapor, and supercritical regions of the phase
    diagram. Implements the equation of state given in:

    R. Tillner-Roth and H. D. Baehr. *An International Standard Formulation for
    The Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for
    Temperatures From 170 K to 455 K and Pressures up to 70 MPa.* J. Phys.
    Chem. Ref. Data, Vol. 23, No. 5, 1994. pp. 657--729.
    http://dx.doi.org/10.1063/1.555958

    For more details, see classes Cantera::PureFluid and tpx::HFC134a in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'HFC-134a')


def CarbonDioxide():
    """
    Create a `PureFluid` object using the equation of state for carbon dioxide.

    The object returned by this method implements an accurate equation of
    state for carbon dioxide that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram. The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances.* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::CarbonDioxide in
    the Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'carbon-dioxide')


def Heptane():
    """
    Create a `PureFluid` object using the equation of state for heptane.

    The object returned by this method implements an accurate equation of
    state for n-heptane that can be used in the liquid, vapor, saturated
    liquid/vapor, and supercritical regions of the phase diagram. The
    equation of state is taken from

    W. C. Reynolds, *Thermodynamic Properties in SI: graphs, tables, and
    computational equations for forty substances.* Stanford: Stanford
    University, 1979. Print.

    For more details, see classes Cantera::PureFluid and tpx::Heptane in the
    Cantera C++ source code documentation.
    """
    return PureFluid('liquidvapor.yaml', 'heptane')
