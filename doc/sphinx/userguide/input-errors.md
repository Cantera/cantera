# Handling Errors in Input Files

During processing of an input file, errors may be encountered. These could be syntax
errors, or could be ones that are flagged as errors by Cantera due to some apparent
inconsistency in the data---an unphysical value, a species that contains an undeclared
element, a reaction that contains an undeclared species, missing species or element
definitions, multiple definitions of elements, species, or reactions, and so on. This
section is intended to help you understand some causes for these errors and their
solutions.

## Syntax Errors

Syntax errors are caught by the YAML parser, and must be corrected before proceeding
further. If a syntax error is encountered, Cantera will raise an exception which
includes the location of the error. Additional information such as a traceback showing
where in the code the input file was being read may be printed as well.

For example, consider the following input file, which is intended to create a gas with
the species and reactions of GRI-Mech 3.0, but is missing the colon which is needed
after the `thermo` key:

```yaml
phases:
- name: gas
  thermo ideal-gas
  kinetics: gas
  elements: [H, O]
  species: [{gri30.yaml/species: all}]
  reactions: [gri30.yaml/reactions]
```

When this definition is imported into an application, an error message like the
following would be printed to the screen, and execution of the program or script would
terminate:

```pytb
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/some/path/cantera/base.pyx", line 25, in cantera._cantera._SolutionBase.__cinit__
    self._init_yaml(infile, phaseid, phases, yaml)
  File "/some/path/cantera/base.pyx", line 49, in cantera._cantera._SolutionBase._init_yaml
    root = AnyMapFromYamlFile(stringify(infile))
cantera._cantera.CanteraError:
***********************************************************************
InputFileError thrown by AnyMap::fromYamlFile:
Error on line 4 of ./gas.yaml:
illegal map value
|  Line |
|     1 | phases:
|     2 | - name: gas
|     3 |   thermo ideal-gas
>     4 >   kinetics: gas
                    ^
|     5 |   elements: [H, O]
|     6 |   species: [{gri30.yaml/species: all}]
|     7 |   reactions: [gri30.yaml/reactions]
***********************************************************************
```

The top part of the error message shows the chain of functions that were called before
the error was encountered. For the most part, these are internal Cantera functions not
of direct concern here. The relevant part of this error message is the part between the
lines of asterisks. This message says that the YAML parser ran into a problem on line 4
of `gas.yaml`. In many cases, including this one, the parser will fail somewhere *after*
the actual problem with the input file, since it must continue parsing until it finds
something that cannot possibly be valid YAML syntax. In this case, the problem from the
parser's perspective is that the key which started on line 3 continues across a new line
before it finds a colon that can be considered as the separator. Since a key can't be
broken across lines like this, the parser indicates the error at the point where it
found the colon. By looking back from the indicated point of the error, we can see that
the problem is the missing colon in the previous line.

## Cantera Errors

Now let's consider the other class of errors, ones that Cantera itself detects.
Continuing the example above, suppose that the missing colon is corrected, and the input
file processed again. Again an error message results, but this time it is from Cantera:

```pytb
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/some/path/cantera/base.pyx", line 25, in cantera._cantera._SolutionBase.__cinit__
    self._init_yaml(infile, phaseid, phases, yaml)
  File "/some/path/cantera/base.pyx", line 49, in cantera._cantera._SolutionBase._init_yaml
    root = AnyMapFromYamlFile(stringify(infile))
cantera._cantera.CanteraError:
***********************************************************************
CanteraError thrown by Phase::addSpecies:
Species 'C' contains an undefined element 'C'.
***********************************************************************
```

The problem is that the phase definition specifies that all species are to be imported
from the `gri30` mechanism, but only the elements H and O are declared. The `gri30`
mechanism contains species composed of the elements H, O, C, N, and Ar. If the
definition is modified to declare these additional elements:

```yaml
phases:
- name: gas
  thermo: ideal-gas
  kinetics: gas
  elements: [H, O, C, N, Ar]
  species: [{gri30.yaml/species: all}]
  reactions: [gri30.yaml/reactions]
```

it may be imported successfully.
