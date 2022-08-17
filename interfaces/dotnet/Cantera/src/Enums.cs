// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Reflection;
using System.Text.RegularExpressions;

namespace Cantera;

/// <summary>
/// Represents a pair of thermodynamic properties that are set or held together.
/// </summary>
public enum ThermoPair
{
#pragma warning disable CS1591 // names are obvious for long constants
    /// <summary>
    /// Density and Pressure
    /// </summary>
    RP, DensityPressure = RP,

    /// <summary>
    /// Temperature and Volume
    /// </summary>
    TV, TemperatureVolume = TV,

    /// <summary>
    /// Enthalpy and Pressure
    /// </summary>
    HP, EnthalpyPressure = HP,

    /// <summary>
    /// Entropy and Pressure
    /// </summary>
    SP, EntropyPressure = SP,

    /// <summary>
    /// Pressure and Volume
    /// </summary>
    PV, PressureVolume = PV,

    /// <summary>
    /// Temperature and Pressure
    /// </summary>
    TP, TemperaturePressure = TP,

    /// <summary>
    /// Internal and Energy
    /// </summary>
    UV, InternalEnergyVolume = UV,

    /// <summary>
    /// Entropy and Temperature
    /// </summary>
    ST, EntropyTemperature = ST,

    /// <summary>
    /// Entropy and Volume
    /// </summary>
    SV, EntropyVolume = SV,

    /// <summary>
    /// Internal and Energy
    /// </summary>
    UP, InternalEnergyPressure = UP,

    /// <summary>
    /// Volume and Enthalpy
    /// </summary>
    VH, VolumeEnthalpy = VH,

    /// <summary>
    /// Temperature and Enthalpy
    /// </summary>
    TH, TemperatureEnthalpy = TH,

    /// <summary>
    /// Entropy and Enthalpy
    /// </summary>
    SH, EntropyEnthalpy = SH,
#pragma warning restore CS1591
}

static class ThermoPairExtensions
{
    static readonly Lazy<IReadOnlyDictionary<ThermoPair, string>> s_interopStringMap =
        new(() => GetThermoPairEnumFieldsWithTwoCharName()
            .ToDictionary(f => (ThermoPair) f.GetValue(null)!, f => f.Name));

    internal static IEnumerable<FieldInfo> GetThermoPairEnumFieldsWithTwoCharName() =>
        typeof(ThermoPair)
            .GetFields(BindingFlags.Static | BindingFlags.Public)
            .Where(f => Regex.IsMatch(f.Name, "^[A-Z]{2}$")); // exactly two uppercase

    public static string ToInteropString(this ThermoPair thermoPair)
    {
        if (s_interopStringMap.Value.TryGetValue(thermoPair, out var interopString))
        {
            return interopString;
        }

        throw new ArgumentOutOfRangeException(nameof(thermoPair));
    }
}

// the constants MUST match what CLIB is expecting

/// <summary>
/// Determines which algorithm is used to find equilibrium.
/// </summary>
public enum EquilibriumSolver
{
    /// <summary>
    /// Allow Cantera to determine the optimum algorithm.
    /// </summary>
    Auto = -1,

    /// <summary>
    /// Solve by using the element potential algorithm.
    /// </summary>
    ElementPotential,

    /// <summary>
    /// Solve by using the general algorithm to minimize Gibbs free energy.
    /// </summary>
    Gibbs,

    /// <summary>
    /// Solved by using the VCS algorithm to minimize Gibbs free energy.
    /// </summary>
    Vcs
}

/// <summary>
/// The
/// </summary>
public enum LogLevel
{
#pragma warning disable CS1591 // names are obvious
    Info,
    Warning,
    Error
#pragma warning restore CS1591
}
