# Writing Examples

Cantera's examples are written in a style so they can be parsed by
[Sphinx Gallery](https://sphinx-gallery.github.io) and rendered nicely on the Cantera
website. This page describes some of the formatting guidelines and helpful syntax for
examples for compatibility with this system.

## Code and reStructuredText Blocks

Each example is divided into a series of blocks that are treated as either code or as
[reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

The first comment block in an example (or the file docstring for Python examples) is
always treated as reStructuredText. It should start with the title of the example
(underlined with `===`), followed by an introduction.

```{literalinclude} ../../../../samples/python/thermo/rankine_units.py
:language: python
:caption: Introduction for a Python example
:start-at: '"""'
:end-before: 'import '
```

```{literalinclude} ../../../../samples/cxx/bvp/blasius.cpp
:language: C++
:caption: Introduction for a C++ example
:start-at: "/*"
:end-at: "*/"
```

Subsequent blocks that should be rendered using reStructuredText syntax are comment
blocks where the first line of the comment (after the comment character) is `%%`. This
special comment block ends at the end of the comment block, or at the first blank line.
A line that is empty except for the comment character indicates a new paragraph within
reStructuredText block.

```{code-block} python
:emphasize-lines: 2-5
# %%
# Plot temperature and major species profiles
# -------------------------------------------
#
# Check and see if all has gone well.

# Extract the spatial profiles as a SolutionArray to simplify plotting specific species
profile = flame.to_array()
```

In the example above, only the highlighted lines are considered as part of the
reStructuredText block. The subsequent comment is treated as a comment in the code
block.

```{literalinclude} ../../../../samples/cxx/bvp/blasius.cpp
:language: C++
:caption: reStructuredText for a C++ example
:start-at: "/* %%"
:end-at: "public:"
```

## Tags / keywords

The initial comment block should also specify a list of tags (keywords) relevant to the
example. The full list of keywords to choose from is defined in
[`example-keywords.txt`](https://github.com/Cantera/cantera/blob/main/doc/example-keywords.txt).

The first tag should be the programming language the example is written in.

New keywords can be added if desired. The purpose of maintaining this list is primarily
to avoid redundant keywords or variant spellings and punctuations.

## Selecting thumbnails

Sphinx Gallery executes the Python examples and renders any figures as part of the
example. By default, the first figure is used as the thumbnail used on the gallery page.
To select a different thumbnail, add a comment in the example using the `sphinx_gallery_thumbnail_number` directive. To pick the last figure, write:

```py
# sphinx_gallery_thumbnail_number = -1
```

## Diagrams

Diagrams can be included to help explain an example. We strongly recommend that diagrams
be drawn as vector artwork (for example, using Inkscape or Illustrator) and saved in the
SVG format to avoid introducing a large number of binary files into our Git repository.

Diagrams should be saved in the directory `doc/sphinx/_static/images/samples`. They
can be included in the HTML version of the example using the following syntax within
a reStructuredText block:

```rst
.. image:: /_static/images/samples/flame-speed.svg
   :width: 50%
   :alt: Freely Propagating Flame
   :align: center
```

## Data files

For examples that require YAML data files besides the existing ones included in the
default Cantera [`data`](https://github.com/Cantera/cantera/tree/main/data) directory,
these files can be added to the
[`cantera-example-data`](https://github.com/cantera/cantera-example-data) repository.
Files in this repository will be installed with Cantera and available under the
`example_data` subdirectory within the standard search path. For example, the mechanism
`n-hexane-NUIG-2015.yaml` can be loaded as:

```py
gas = ct.Solution("example_data/n-hexane-NUIG-2015.yaml")
```

When developing examples, the use of well-documented mechanisms that have been published
in the peer-reviewed literature is highly encouraged.

## Other Style Guidelines

- Use the following markers to underline headings:
  - Example title with `====`
  - Section headings with `---`
  - Subsections with `^^^`
