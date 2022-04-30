# Installing Cantera

## Basic Instructions

To compile using the default options, run `scons build` followed by
`scons install`.

Configuration options are specified with `name=value` on the command line, for example:
`scons build optimize=n prefix=/home/$USER/cantera`

The full list of configuration options and their default values can be shown by running
`scons help --options`. The list of available `scons` commands (such as `build`) can be
shown by running `scons help`.


## Detailed Instructions

See the instructions available [online](https://cantera.org/install/compiling-install.html).
