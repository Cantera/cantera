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
    readonly static Lazy<IReadOnlyDictionary<ThermoPair, String>> InteropStringMap =
        new(() => typeof(ThermoPair)
            .GetFields(BindingFlags.Static | BindingFlags.Public)
            .Where(f => Regex.IsMatch(f.Name, "^[A-Z]{2}$")) // match exactly two uppercase
            .ToDictionary(f => (ThermoPair) f.GetValue(null)!, f => f.Name));
    
    public static string ToInteropString(this ThermoPair thermoPair)
    {
        if(InteropStringMap.Value.TryGetValue(thermoPair, out var interopString))
            return interopString;

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
