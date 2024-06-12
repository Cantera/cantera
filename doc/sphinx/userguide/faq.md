# Frequently Asked Questions

```{hint}
If your question isn’t answered here, consider asking us on the
[Cantera Users' Group](https://groups.google.com/g/cantera-users).
```

(sec-faq-installation)=
## Installation

:::{dropdown} How do I determine the cause of the error "DLL load failed while importing _cantera: The specified module can not be found."?
:color: primary

1. Install the [Dependencies](https://github.com/lucasg/Dependencies) tool.
2. *If you are using Conda*: Run this tool from a terminal where your Cantera Conda
   environment is active. That is, run the command:
   ```
   start "C:\Program Files\Dependencies\DependenciesGUI.exe"
   ```
   from that terminal, replacing the path with the path to `DependenciesGUI.exe` on your
   computer. This step is necessary to have the `PATH` environment variable set
   correctly.
3. Use Dependencies to open the Cantera Python extension module. If you installed
   Mambaforge at `C:\mambaforge` and Cantera in an environment named `ct-env`, then this
   file would be named something like
   `C:\mambaforge\envs\ct-env\Lib\site-packages\cantera\_cantera.cp312-win_amd64.pyd`.
4. Identify any DLLs that the tool is unable to find, or if there are methods that are
   not being resolved. One thing to keep an eye out for is if there are any libraries
   that are being loaded from unexpected locations (for example, other than from your
   conda installation and the Windows system directories).
:::

## Input Files

:::{dropdown} Where can I find chemical mechanisms to use with Cantera?
:color: primary

There are a few sites that distribute mechanisms in the Cantera YAML format:
- The [Caltech Explosion Dynamics Laboratory](https://shepherd.caltech.edu/EDL/PublicResources/sdt/cti_mech.html)
  provides a number of mechanisms for combustion applications.
- [CollectionOfMechanisms](https://github.com/jiweiqi/CollectionOfMechanisms) is a
  user-maintained GitHub repository of mechanisms that have been obtained from
  scientific publications and other sources.

Many research groups maintain pages with mechanisms they developed, provided in the
  Chemkin format. These mechanisms can be converted to the Cantera YAML format using the
  [`ck2yaml`](/userguide/ck2yaml-tutorial) tool. The following is an incomplete list:

- The [University of Galway Combustion Chemistry Centre](https://www.universityofgalway.ie/combustionchemistrycentre/mechanismdownloads/)
  provides mechanisms for combustion
- [DETCHEM](https://www.detchem.com/mechanisms) provides a variety of mechanisms for
  catalytic systems that were developed by the research group of Olaf Deutschmann at the
  Karlsruhe Institute of Technology.
- [Lawrence Livermore National Laboratory ](https://combustion.llnl.gov/mechanisms)
  provides a number of combustion mechanisms

```{caution}
The inclusion of a site in the above lists does not constitute an endorsement of the
chemical mechanisms provided there. You must use your own judgement to determine if a
mechanism is appropriate for any particular scientific or engineering purpose.
```

:::


:::{dropdown} How do I fix errors converting Chemkin input files to Cantera's format?
:color: primary

See our documentation on [](sec-debugging-chemkin). If you're encountering an issue not
described there, please post a question on the
[Cantera Users' Group](https://groups.google.com/g/cantera-users).
:::

## Thermodynamics & Equilibrium

:::{dropdown} Why don't the enthalpy and internal energy values calculated by Cantera match data from other sources (CoolProp / NIST / steam tables)? Why does my mixture have a negative enthalpy?
:color: primary

Reference states in thermodynamics are arbitrary, since the observable quantities of
interest (work done, heat transferred) do not depend on absolute quantities but on
differences between states. Different data sources use different conventions for the
reference state. Cantera uses the convention commonly used in reaction thermodynamics,
where pure elements in their natural state (for example, oxygen as O₂ gas and carbon as
graphite) have zero enthalpy at a temperature of 298.15 K and a pressure of 1 atm. This
convention provides a consistent treatment for the change in enthalpy associated with
chemical reactions. In contrast, it is common for tabulated data for a single pure
substance to use a different reference state, such as CO₂ where the enthalpy may be
referenced to the saturated liquid at 273.15 K.
:::

:::{dropdown} Why doesn't the heat capacity calculated by Cantera match what I get from CEA?
:color: primary

Cantera’s definition of $c_p$ is the partial derivative of enthalpy with respect to
temperature, with pressure and composition held constant. The value returned by CEA is
the derivative while holding the system at chemical equilibrium. The former definition
is the one that appears in governing equations for multi-species systems, such as
[well-stirred reactors](/reference/reactors/ideal-gas-constant-pressure-reactor). While
the latter property is of interest in certain physical situations, it is of somewhat
limited use computationally, since it is *only* defined for mixtures at equilibrium.

You can use Cantera to compute a value like what’s returned by CEA by computing $dh/dT$
using a finite difference method where the mixture is equilibrated at constant $T$ and
$P$ at each point before calculating the enthalpy.
[`sound_speed.py`](/examples/python/thermo/sound_speed) presents a related example.

CEA provides a value consistent with Cantera's definition which it calls the $c_p$ "with
frozen reactions."
:::

:::{dropdown} Why doesn't changing reactions affect the equilibrium state calculated by Cantera? How does Cantera find the equilibrium state even when no reactions are defined?
:color: primary
Chemical equilibrium is not defined in terms of reactions. For whatever set of species
are allowed in a phase or multi-phase mixture, it is defined as the composition that
minimizes a particular thermodynamic state variable while holding two others constant
(for example, for equilibrium at constant temperature and pressure, the Gibbs free
energy is minimized). This optimization is independent of any reactions that are
defined, and the resulting composition may be different from what you would get, for
example, from integrating a reactor network using a set of irreversible chemical
reactions or a set of reactions that does not form a complete basis set for the defined
species.
:::

## Reactor Networks

:::{dropdown} How do I set the residence time for a reactor?
:color: primary

Cantera defines reactor flows in terms of inlet and outlet mass flow rates, with
residence time being an output that can be calculated:

$$ t_\t{res} = \frac{m_\t{combustor}}{\dot{m}_\t{in}} $$

To achieve a specified residence time, you can define the mass flow rate to be a
function of the combustor state. For example, the following approach is used in the
example [`combustor.py`](/examples/python/reactors/combustor):

```py
def mdot(t):
    return combustor.mass / residence_time

inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)
```
:::

:::{dropdown} How do I understand errors from ReactorNet that give a list of "components with largest weighted error estimates"?
:color: primary

Calls to `ReactorNet.step` or `ReactorNet.advance` may result in error messages like the
following:
```
 CanteraError thrown by CVodesIntegrator::integrate:
CVodes error encountered. Error code: -4
At t = 0.151127 and h = 1.61901e-09, the corrector convergence test failed repeatedly or with |h| = hmin.

Components with largest weighted error estimates:
834: 297.28141792959997
762: -273.0999483219796
8: 1.8673611455956203
99: -0.018812031839692638
100: 0.01850247743429419
2: -0.002587710808893717
771: 0.0007865229822411536
843: -0.0007861126902119126
184: -0.0006244451494927259
231: 0.0006129536140016525
```

The list labeled "Components with largest weighted error estimates" provides a list of
component indices followed by their corresponding error estimates. Large error estimates
tend to be associated with rapidly changing variables which require very small time
steps. If you are using the Python interface and the `ReactorNet` object is named `net`,
you can identify the variables associated with these indices:

```pycon
>>> net.component_name(834)
'ABC'
>>> net.component_name(762)
'XYZ'
```

It is sometimes the case that a reaction mechanism will contain one or more reactions
where at some temperatures, either the forward or reverse rate constant becomes high and
nonphysical, particularly if the mechanism has not been designed for use at in a
particular temperature range. For this example, you can find reactions involving the
species with the highest errors, assuming a `Solution` object named `gas`:

```pycon
>>> for i, R in enumerate(gas.reactions()):
...     spec = R.reactants | R.products
...     if 'ABC' in spec and 'XYZ' in spec:
...         print(i, R)
244 ABC <=> XYZ
```

For this example, resolution of the problem would then involve investigating the rate
parameterization for reaction 244 or the thermodynamic data for species ABC and XYZ,
which determine the reverse rate constant.

In other cases, the time reached by the integrator (here, `t = 0.151127`) can provide a
hint at the source of the problem. For example, if this time is near the time of a
discontinuity in the inputs to a reactor network, such as a valve opening or closing,
then the problem might be solved by using a smoother time function or decreasing the
maximum integrator timestep.
:::

## 1D Reacting Flows


:::{dropdown} Why can't I calculate the flame speed for a mixture with an inlet temperature of ~1000 K or higher?
:color: primary

Two failure modes are common under these conditions:

1. when the `auto=True` option is specified to `FreeFlame.solve`, the solver keeps
expanding the domain and never finds a solution satisfying its tolerances
2. the solver converges but the solution is a strong function of the distance from the
inlet to the flame.

The cause of both these problems is that at high inlet temperatures, reactions are
already starting to occur at that inlet temperature. This violates the standard boundary
conditions for the laminar flame problem, where the reaction rates need to go to zero
for the inlet mixture. Here, however, finite rates at an inlet that is a finite distance
from the flame mean that the mixture has changed before it even reaches the flame, and
the computed flame speed becomes a function of that distance. Cantera's solver (with the
"auto" option enabled) tries to keep the inlet boundary far enough away from the flame
to avoid problems with non-zero diffusive fluxes across the inlet, and that leads it to
further widening the domain as the temperature continues to go up.

Considering a domain with the temperature fixed point at a distance $d$ from the inlet,
the ignition delay time $\tau_\t{ig}$ sets a lower bound on the calculated flame speed.
That is, even in the absence of any diffusion of heat or radicals into the unburned
mixture, the calculated flame speed will be $S_u \approx d / \tau_\t{ig}$ or, accounting
for transport, $S_u > d / \tau_\t{ig}$. A correct flame speed calculation can only be
obtained when $S_u$ is independent of $d$ and $d$ is larger than the flame thickness.
:::
