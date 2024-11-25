(sec-yaml-elements)=
# Elements

`element` entries are needed only when defining custom elements that are not standard
chemical elements, or defining specific isotopes.

The fields of an `element` entry are:

`symbol`
: The symbol used for the element, as used when specifying the composition of species.

`atomic-weight`
: The atomic weight of the element, in unified atomic mass units (dalton).

`atomic-number`
: The atomic number of the element. Optional.

`entropy298`
: The standard molar entropy of the element at 298.15 K. Optional.

An example `elements` section:

```yaml
elements:
- symbol: C13
  atomic-weight: 13.003354826
  atomic-number: 6
- symbol: O-18
  atomic-weight: 17.9991603
```
