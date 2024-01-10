(sec-thermobuild)=
# Converting Data from NASA ThermoBuild

% TODO: Link to NASA9 Science docs

Thermodynamic data for a range of species can be obtained from the
[NASA ThermoBuild](https://cearun.grc.nasa.gov/ThermoBuild/index_ds.html) tool.
This thermodynamic data is described using the 9-coefficient NASA polynomial parameterization, which is implemented by Cantera.

To generate an input file containing thermodynamic data, use the ThermoBuild web
interface to select the elements you want to include, and click *process*. Then, select
the species you want to include in your mechanism and click *continue*. Select the
content on the resulting page starting with the line that says `thermo`, and ending with
the line that says `END REACTANTS` and save it to a file.

Next, we need to make a few modifications to this file so it contains all of the
information needed to be processed by Cantera's {term}`CK` file parser. For this
example, let's suppose you have selected just the species `CO` and `CO2` from
ThermoBuild.

- Modify the first line of the file so that it reads `thermo nasa9`. This is done to
  distinguish this input format from the one used for the 7-coefficient NASA
  polynomials.
- Replace the last two lines (`END PRODUCTS` and `END REACTANTS`) with a single line
  that reads `END`.

## Creating species definitions only

If you only want to generate a species database that can be referenced from other
Cantera input files, the above modifications are sufficient. The last step is simply to
convert the input file to the Cantera YAML format using `ck2yaml`. If you named the file
`mythermo.txt`, then the command to convert it would be:

```bash
ck2yaml --thermo=mythermo.txt
```

This will generate the file `mythermo.yaml`. You can then reference the species
definitions from this file in phase definitions in other Cantera input files, for
example:

```yaml
phases:
- name: gas
  thermo: ideal-gas
  species:
  - {mythermo.yaml/species: all}
```

## Creating a complete phase definition

To create a complete phase definition, you also need to add two sections to the top of
the ThermoBuild input file. First, a section declaring all of the elements:

```
elements
C O
end
```

And second, a section declaring all of the species:

```
species
CO CO2
end
```

The resulting input file should look like the following:

```
elements
C O
end

species
CO CO2
end

thermo nasa9
   200.000  1000.000  6000.000 20000.000   9/09/04
CO                Gurvich,1979 pt1 p25 pt2 p29.
 3 tpis79 C   1.00O   1.00    0.00    0.00    0.00 0   28.0101000    -110535.196
    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
 1.489045326D+04-2.922285939D+02 5.724527170D+00-8.176235030D-03 1.456903469D-05
-1.087746302D-08 3.027941827D-12                -1.303131878D+04-7.859241350D+00
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
 4.619197250D+05-1.944704863D+03 5.916714180D+00-5.664282830D-04 1.398814540D-07
-1.787680361D-11 9.620935570D-16                -2.466261084D+03-1.387413108D+01
   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
 8.868662960D+08-7.500377840D+05 2.495474979D+02-3.956351100D-02 3.297772080D-06
-1.318409933D-10 1.998937948D-15                 5.701421130D+06-2.060704786D+03
CO2               Gurvich,1991 pt1 p27 pt2 p24.
 3 g 9/99 C   1.00O   2.00    0.00    0.00    0.00 0   44.0095000    -393510.000
    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
 4.943650540D+04-6.264116010D+02 5.301725240D+00 2.503813816D-03-2.127308728D-07
-7.689988780D-10 2.849677801D-13                -4.528198460D+04-7.048279440D+00
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
 1.176962419D+05-1.788791477D+03 8.291523190D+00-9.223156780D-05 4.863676880D-09
-1.891053312D-12 6.330036590D-16                -3.908350590D+04-2.652669281D+01
   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
-1.544423287D+09 1.016847056D+06-2.561405230D+02 3.369401080D-02-2.181184337D-06
 6.991420840D-11-8.842351500D-16                -8.043214510D+06 2.254177493D+03
END
```

This file (saved for example as `myphase.txt`) can then be converted to the Cantera YAML
format using the [`ck2yaml`](/yaml/ck2yaml) utility from a shell:

```bash
ck2yaml --input=myphase.txt
```

This will generate a an input file named `myphase.yaml` with a phase named `gas` that
can be directly imported in Cantera.
