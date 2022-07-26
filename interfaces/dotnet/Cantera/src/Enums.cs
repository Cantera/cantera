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
#pragma warning disable CS1591 // name are obvious
    RP, DensityPressure = RP,

    TV, TemperatureVolume = TV,

    HP, EnthalpyPressure = HP,

    SP, EntropyPressure = SP,

    PV, PressureVolume = PV,

    TP, TemperaturePressure = TP,

    UV, InternalEnergyVolume = UV,

    ST, EntropyTemperature = ST,

    SV, EntropyVolume = SV,

    UP, InternalEnergyPressure = UP,

    VH, VolumeEnthalpy = VH,

    TH, TemperatureEnthalpy = TH,

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
/// Determines which algorithm is used to find equilibirum.
/// </summary>
public enum EquilibriumSolver
{
    /// <summary>
    /// Allow Cantera to determine the optimum algorithm.
    /// </summary>
    Auto = -1,

    /// <summary>
    /// Solve by using elment potential algorithm.
    /// </summary>
    ElementPotential,

    /// <summary>
    /// Solve by using the general alogrithm to minimize Gibbs free energy.
    /// </summary>
    Gibbs,

    /// <summary>
    /// Solved by using the VCS alogrithm to minimize Gibbs free energy.
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
