// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Diagnostics.CodeAnalysis;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Represents a thermodynamic phase.
/// </summary>
public partial class ThermoPhase
{
    [SuppressMessage("Usage", "CA2213: Disposable field not disposed",
        Justification = "Field actually is disposed, in ExtraDispose()")]
    readonly SolutionHandle _sol;

    readonly Lazy<SpeciesCollection> _species;

    /// <summary>
    /// The collection of species that make up this phase.
    /// </summary>
    public SpeciesCollection Species => _species.Value;

    internal ThermoPhase(string filename, string? phaseName)
    {
        _sol = LibCantera.sol_newSolution(filename, phaseName ?? "", "none");
        _handle = LibCantera.sol_thermo(_sol);

        _species = new(() => new SpeciesCollection(_handle));
    }

    /// <summary>
    /// Simulates bringing the phase to thermodynamic equilibrium by holding the
    /// specified <see cref="ThermoPair" /> constant and using the algorithm(s)
    /// identified by the solver string.
    /// </summary>
    public void Equilibrate(ThermoPair thermoPair,
                            string solver = "auto",
                            double tolerance = 1e-9, int maxSteps = 1000,
                            int maxIterations = 100, int logVerbosity = 0)
    {
        var interopString = thermoPair.ToInteropString();

        LibCantera.thermo_equilibrate(_handle, interopString, solver,
            tolerance, maxSteps, maxIterations, logVerbosity);
    }

    /// <summary>
    /// Sets the given pair of thermodynamic properties for this phase together.
    /// </summary>
    public void SetPair(ThermoPair pair, double first, double second) =>
        _ = pair switch
        {
            ThermoPair.DP => LibCantera.thermo_setState_DP(_handle, first, second),
            ThermoPair.TV => LibCantera.thermo_setState_TV(_handle, first, second),
            ThermoPair.HP => LibCantera.thermo_setState_HP(_handle, first, second),
            ThermoPair.SP => LibCantera.thermo_setState_SP(_handle, first, second),
            ThermoPair.PV => LibCantera.thermo_setState_PV(_handle, first, second),
            ThermoPair.TP => LibCantera.thermo_setState_TP(_handle, first, second),
            ThermoPair.UV => LibCantera.thermo_setState_UV(_handle, first, second),
            ThermoPair.ST => LibCantera.thermo_setState_ST(_handle, first, second),
            ThermoPair.SV => LibCantera.thermo_setState_SV(_handle, first, second),
            ThermoPair.UP => LibCantera.thermo_setState_UP(_handle, first, second),
            ThermoPair.VH => LibCantera.thermo_setState_VH(_handle, first, second),
            ThermoPair.TH => LibCantera.thermo_setState_TH(_handle, first, second),
            ThermoPair.SH => LibCantera.thermo_setState_SH(_handle, first, second),
            _ => throw new ArgumentOutOfRangeException(nameof(pair))
        };

    partial void ExtraDispose()
    {
        _sol.Dispose();
    }
}
