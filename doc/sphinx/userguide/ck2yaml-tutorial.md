# Converting Chemkin Format Files

Many existing reaction mechanism files are in **{term}`CK` format**, by which we mean
the input file format developed for use with the Chemkin-II software package (and
subsequent releases) as specified in the report describing the Chemkin software
{cite:p}`kee1989`. Cantera comes with a converter utility program `ck2yaml` (or
`ck2yaml.py`) that converts CK format into Cantera's YAML format. If you want to convert
a Chemkin-format file to YAML format, or you're having errors when you try to do so,
this section will help.

## Converting gas-phase mechanisms

:::{seealso}
For documentation of all the command line options that can be used with `ck2yaml`, see
its [documentation in the reference section](sec-ck2yaml).
:::

To convert a gas phase mechanism where the thermodynamic and transport data are provided
as separate files, run the following command from a terminal (command prompt):

```bash
ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat
```

If the `ck2yaml` script is not on your path but the Cantera Python module is,
`ck2yaml` can also be used by running:

```bash
python -m cantera.ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat
```

:::{tip}
If you're using Cantera from IPython or a Jupyter Notebook, you can run the above shell
command from Python by prefixing it with an exclamation point. For example:

```ipython
!python -m cantera.ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat
```
:::

## Converting surface mechanism files

For mechanisms involving both surface and gas phase species, an additional input file
is needed to define the surface species and reactions occurring on the surface. This
file is specified with the `--surface` option:

```bash
python -m cantera.ck2yaml --input=gas.inp --thermo=therm.dat --surface=surf.inp
```

## Converting standalone thermo data

An input file containing only species definitions (which can be referenced from
phase definitions in other input files) can be created by specifying only a
thermo file. For example:

```bash
ck2yaml --thermo=therm.dat
```

## Debugging common errors in CK files

When `ck2yaml` encounters an error, it attempts to print the surrounding information to
help you to locate the error. Many of the most common errors are due to an inconsistency
of the input files from their standard, as defined in the report for Chemkin referenced
above. Here, we describe some requirements of the CK format that are common sources of
errors.

:::{tip}
Many existing CK format files cause errors in `ck2yaml` when they are processed. Some of
these errors may be avoided by specifying the `--permissive` option. This option allows
certain recoverable parsing errors (for example, duplicate transport or thermodynamic
data) to be ignored.
:::

### Input file sections

Each section of a CK input file must start with a keyword representing that section and
end with the keyword `END`. Keywords that may begin a section include:

  - `ELEMENTS` or `ELEM`
  - `SPECIES` or `SPEC`
  - `THERMO` or `THERMO ALL`
  - `REACTIONS` or `REAC`
  - `TRANSPORT`

### Thermo data

The thermodynamic data is read in a fixed format. This means that each column of the
input has a particular meaning. The format for each thermodynamic entry should be as
follows:

```
CH3OH             L 8/88C   1H   4O   1     G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
```

:::{caution}
Many common errors are generated because information is missing or in the wrong column.
Check thoroughly for extraneous or missing spaces. Also make sure that the input file
does not contain tab characters.
:::

The following table is adapted from the Chemkin manual {cite:p}`kee1989` to describe
the column positioning of each required part of the entry. Empty columns should be
filled with spaces.

:::{list-table}
:header-rows: 1

* - Line No.
  - Columns
  - Contents
* - 1
  - 1--18
  - Species Name
* - 1
  - 19--24
  - Date (Optional)
* - 1
  - 25--44
  - Atomic Symbols and formula
* - 1
  - 45
  - Phase of species (`S`, `L`, or `G`)
* - 1
  - 46--55
  - Low temperature
* - 1
  - 56--65
  - High temperature
* - 1
  - 66--73
  - Common temperature
* - 1
  - 74--78
  - Additional Atomic Symbols
* - 1
  - 80
  - The integer `1`
* - 2
  - 1--75
  - Coefficients $a_1$ to $a_5$ for the upper temperature interval
* - 2
  - 80
  - The integer `2`
* - 3
  - 1--75
  - Coefficients $a_6,\ a_7$ for the upper temperature interval, and $a_1,\ a_2,\ a_3$
    for the lower temperature interval
* - 3
  - 80
  - The integer `3`
* - 4
  - 1--60
  - Coefficients $a_4$ through $a_7$ for the lower temperature interval
* - 4
  - 80
  - The integer `4`
:::

#### Line 1

The first 18 columns are reserved for the species name. The name assigned to the
species in the thermodynamic data must be the same as the species name defined in the
`SPECIES` section. If the species name is shorter than 18 characters, the rest of the
characters should be filled by spaces.

The next six columns (columns 19--24) are typically used to write a date; they are not
used further.

The next 20 columns (25--44) are used to specify the elemental composition of the
species. Five characters are used for each element: the first two specify the element
symbol, while the remaining three specify the number of atoms of that element in the
species. Up to four elements may be specified this way. For species with more than four
elements, the [extended composition format](sec-ck-extended-composition) described below
may be used.

In column 45, the phase of the species (`S`, `L`, or `G` for solid, liquid, or gas
respectively) should be specified.

The next 28 columns are reserved for the temperatures that delimit the ranges of the
polynomials specified on the next several lines. The first two temperatures have a width
of 10 columns each (46--55 and 56--65), and represent the lowest temperature and highest
temperature for which the polynomials are valid. The last temperature has a width of 8
columns (66--73) and is the **common** temperature, where the switch from low to high
occurs.

The next 5 columns (74--78) are reserved for atomic symbols and are usually left blank
for the default behavior.

Column 79 is blank; and finally, the row is ended in column 80 with the integer `1`.

#### Lines 2--4

The next three lines of the thermodynamic entry have a similar format. They contain the
coefficients of the [7-coefficient polynomial formulation](sec-thermo-nasa7) in two
temperature regions.

The second row of the thermo entry (the first after the information row) contains the
first five coefficients that apply to the temperature range between the midpoint and the
upper limit. 15 columns are alloted for each coefficient (for a total of 75 columns),
with no spaces between them. Although the entry above shows spaces between positive
coefficients, it is to be noted that this is done only for formatting consistency with
other lines that contain negative numbers. After the coefficients, four spaces in
columns 76--79 are followed by the integer `2` in column 80.

On the next line, the last two coefficients for the upper temperature range and the
first three coefficients for the lower temperature range are specified. Once again, this
takes up the first 75 columns, columns 76--79 are blank, and the integer `3` is in
column 80.

Finally, on the last line of a particular entry, the last four coefficients of the lower
temperature range are specified in columns 1--60, 19 blank spaces are present, and the
integer `4` is in column 80.

(sec-ck-extended-composition)=
#### Extended composition format

If the number of atoms of an element in a thermodynamic entry has more than 3 digits,
the standard format for the composition cannot be used. To allow for such species, the
following extended format is allowed. The element symbol should have a `0` in the first
line of the entry. An ampersand (`&`) is added after the index of the first line, and
the element symbols and their amounts should be written on the next line as follows:

```
BIN6J      PYRENEJ1     C   0H   0    0    0G   300.000  5000.000 1401.000     1&
C 778    H 263
 3.63345177E+01 3.13968020E-02-1.09044660E-05 1.71125597E-09-1.00056355E-13    2
 4.05143093E+04-1.77494305E+02-1.20603441E+01 1.59247554E-01-1.41562602E-04    3
 6.26071650E-08-1.09305161E-11 5.56473533E+04 7.68451211E+01                   4
```

or on separate lines with ampersand (`&`) as the last character on the line:

```
BIN6       PYRENE       C   0H   0    0    0G   300.000  5000.000 1401.000     1&
C      778&
H      264
 3.65839677E+01 3.36764102E-02-1.16783938E-05 1.83077466E-09-1.06963777E-13    2
 9.29809483E+03-1.81272070E+02-1.29758980E+01 1.63790064E-01-1.43851166E-04    3
 6.31057915E-08-1.09568047E-11 2.48866399E+04 7.94950474E+01                   4
```

### Transport parameters

The transport data file also has a specified format, as described in {cite:p}`kee1986`,
although the format is not as strict as for the thermodynamic entries. In particular,
the first 15 columns of a line are reserved for the species name. *One common source of
errors is a species that is present in the transport data file, but not in the
thermodynamic data or in the species list; or a species that is present in the species
list but not the transport data file.* The rest of the columns on a given line have no
particular format, but must be present in the following order:

:::{list-table}
:header-rows: 1

* - Parameter Number
  - Parameter Name
* - 1
  - An integer with value 0, 1, or 2 indicating monatomic, linear, or non-linear
    molecular geometry.
* - 2
  - The Lennard-Jones potential well depth $\varepsilon/k_B$ in Kelvin
* - 3
  - The Lennard-Jones collision diameter $\sigma$ in Angstrom
* - 4
  - The dipole moment $\mu$ in Debye
* - 5
  - The polarizability $\alpha$ in Angstroms cubed
* - 6
  - The rotational relaxation collision number $Z_{rot}$ at 298 K
:::

Another common error is if all six of these numbers are not present for every species.

### Other formatting issues

It may be the case that scientific formatted numbers are missing the `E`. In this case,
numbers often show up as `1.1+01`, when they should be `1.1E+01`. You can fix this with
a Regular Expression "find and replace":

```
Find: (\d+\.\d+)([+-]\d+)
Replace: $1E$2
```
