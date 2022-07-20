// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Reflection;
using System.Text.RegularExpressions;

namespace Cantera;

public enum ThermoPair
{
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
}

static class ThermoPairExtensions
{
    internal static IEnumerable<FieldInfo> GetThermoPairEnumFieldsWithTwoCharName() =>
        typeof(ThermoPair)
            .GetFields(BindingFlags.Static | BindingFlags.Public)
            .Where(f => Regex.IsMatch(f.Name, "^[A-Z]{2}$")); // exactly two uppercase

    readonly static Lazy<IReadOnlyDictionary<ThermoPair, string>> InteropStringMap =
        new(() => GetThermoPairEnumFieldsWithTwoCharName()
            .ToDictionary(f => (ThermoPair) f.GetValue(null)!, f => f.Name));

    public static string ToInteropString(this ThermoPair thermoPair)
    {
        if (InteropStringMap.Value.TryGetValue(thermoPair, out var interopString))
        {
            return interopString;
        }

        throw new ArgumentOutOfRangeException(nameof(thermoPair));
    }
}

// the constants MUST match what CLIB is expecting
public enum EquilibriumSolver
{
    Auto = -1,

    ElementPotential,

    Gibbs,

    Vcs
}

public enum LogLevel
{
    Info,
    Warning,
    Error
}
