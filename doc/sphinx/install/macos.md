(sec-install-macos)=
# macOS Matlab Toolbox

The Cantera Matlab toolbox can be installed using a macOS-specific installer. The
toolbox requires macOS version 10.15 (Catalina) or higher and a copy of Matlab built for
Intel processors (which will run under Rosetta 2 on ARM-based processors).

:::{seealso}
To install the Cantera Python package, see the [pip](pip) or [conda](conda)
instructions. The Python package is required if:

- You need to convert legacy input files to YAML
- You need to convert Chemkin-format input files to YAML
:::

## Installation

Download the Matlab Interface Installer package from
[GitHub](https://github.com/Cantera/cantera/releases/)

When the file has downloaded, find it in Finder, hold *Control* and click the file.
Choose *Open* from the resulting menu, and select *Open* in the security dialog that
appears. Click *Continue* to proceed in the installer (noting that the installer may
open in the background; you can find its icon on the Dock), agreeing to the
[Cantera license terms](https://github.com/Cantera/cantera/blob/main/License.txt) and
the terms of the other open source software that we use.

By default, the installer will add some lines to the file
`$HOME/Documents/MATLAB/startup.m` to enable loading the Cantera toolbox when Matlab
starts. If you wish to disable this, click *Customize* and de-select the *Install
startup.m script* option. Finally, clicking *Install* will install the interface to the
`$HOME/Applications/Cantera` folder.

## Testing the installation

Open Matlab and enter the following code:

```matlab
gas = Solution('gri30.yaml')
h2o = Solution('liquidvapor.yaml','water')
```
