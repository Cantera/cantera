// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
    DP, DensityPressure = DP,

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
    public static string ToInteropString(this ThermoPair thermoPair) =>
        thermoPair switch
        {
            ThermoPair.DP => "DP",
            ThermoPair.TV => "TV",
            ThermoPair.HP => "HP",
            ThermoPair.SP => "SP",
            ThermoPair.PV => "PV",
            ThermoPair.TP => "TP",
            ThermoPair.UV => "UV",
            ThermoPair.ST => "ST",
            ThermoPair.SV => "SV",
            ThermoPair.UP => "UP",
            ThermoPair.VH => "VH",
            ThermoPair.TH => "TH",
            ThermoPair.SH => "SH",
            _ => throw new ArgumentOutOfRangeException(nameof(thermoPair))
        };
}

// the constants MUST match what CLib is expecting

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
