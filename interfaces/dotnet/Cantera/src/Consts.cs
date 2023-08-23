// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera;

/// <summary>
/// Contains well known thermodynamic constants.
/// </summary>
public static class Consts
{
    /// <summary>
    /// Avogadro's Number @f$ N_{\mathrm{A}} @f$ [number/kmol]
    /// </summary>
    public const double Avogadro = 6.02214076e26;

    /// <summary>
    /// Boltzmann constant @f$ k @f$ [J/K]
    /// </summary>
    public const double Boltzmann = 1.380649e-23;

    /// <summary>
    /// Planck constant @f$ h @f$ [J-s]
    /// </summary>
    public const double Planck = 6.62607015e-34;

    /// <summary>
    /// Elementary charge @f$ e @f$ [C]
    /// </summary>
    public const double ElectronCharge = 1.602176634e-19;

    /// <summary>
    /// Speed of Light in a vacuum @f$ c @f$ [m/s]
    /// </summary>
    public const double LightSpeed = 299792458.0;

    /// <summary>
    /// Electron Mass @f$ m_e @f$ [kg]
    /// </summary>
    public const double ElectronMass = 9.1093837015e-31;

    /// <summary>
    /// Universal Gas Constant @f$ R_u @f$ [J/kmol/K]
    /// </summary>
    public const double GasConstant = Avogadro * Boltzmann;

    /// <summary>
    /// Faraday constant @f$ F @f$ [C/kmol]
    /// </summary>
    public const double Faraday = ElectronCharge * Avogadro;

    /// <summary>
    /// Stefan-Boltzmann constant @f$ \sigma @f$ [W/m2/K4]
    /// </summary>
    public const double StefanBoltzmann = 5.670374419e-8;

    /// <summary>
    /// One atmosphere [Pa]
    /// </summary>
    public const double OneAtm = 1.01325e5;

    /// <summary>
    /// One bar [Pa]
    /// </summary>
    public const double OneBar = 1.0E5;
}
