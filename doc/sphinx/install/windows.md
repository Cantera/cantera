(sec-install-windows)=
# Windows Packages

Windows installers are provided for stable versions of Cantera. These installers
provide the Matlab toolbox and header/library files that can be used to compile
C++ applications.

:::{seealso}
To install the Cantera Python package, see the [pip](pip) or [conda](conda)
installation instructions. The Python package is required if:

- You need to convert Chemkin-format input files to YAML
- You need to convert legacy CTI or XML input files to YAML
:::

1. **Remove old versions of Cantera**

- Use The Windows "Add/Remove Programs" interface to remove previous versions of
  the `Cantera` package.

2. **Install Cantera**

- Go to the [Cantera Releases](https://github.com/Cantera/cantera/releases) page and
  download **Cantera-3.0.0-x64.msi**.
- Run the installer and follow the prompts.

3. **Configure Matlab**

- Launch Matlab
- Go to *File->Set Path...*
- Select *Add with Subfolders*
- Browse to the folder `C:\Program Files\Cantera\matlab\toolbox`
- Select *Save*, then *Close*.

4. **Test the installation**

- From the Matlab prompt, run:

  ```matlab
  gas = Solution('gri30.yaml')
  h2o = Solution('liquidvapor.yaml', 'water')
  ```
