# Making a New Cantera Release

Checklist for making a new Cantera release. This is mainly useful to Cantera
maintainers.

## Before the Release

- Reserve DOI on Zenodo.

- Check the [documentation](https://cantera.org/dev/) for accessibility concerns, for
  example by using the
  [Lighthouse](https://developer.chrome.com/docs/lighthouse/overview) tool built into
  Chrome.
  - Check both Sphinx-generated and Doxygen-generated pages.
  - In addition to the automated checks done by Lighthouse, consider some of the manual
    checks that it suggests, such as testing that pressing the "tab" key cycles through
    page elements in a logical order.
  - The most likely changes affecting accessibility would come from changes to the
    versions of the tools used to build the documentation, particularly Doxygen, Sphinx,
    and the Pydata Sphinx theme. For consistency, these all have pinned versions used
    in the `docs` CI build defined in `.github/workflows/main.yml`.

- Check the versions of dependencies that we suggest work for installing/compiling
  Cantera.
  - Installation instructions in `doc/sphinx/install/*.md`
  - `doc/sphinx/develop/compiling/dependencies.md`
  - `doc/sphinx/develop/compiling/compilation-reqs.md`

## Release Pull Request

- Make sure version is set correctly:
  - `SConstruct`: `'cantera_version'` (most widely used)
  - `README.rst`
  - `doc/doxygen/Doxyfile`: `'PROJECT_NUMBER'`
  - `ck2yaml.py`, `cti2yaml.py`, and `ctml2yaml.py`
  - For final (non alpha/beta) releases, make sure any alpha/beta suffixes are removed
    from the `cantera-version` field of any YAML input files that have been updated
  - Installation instructions (`doc/sphinx/install/*.md`)

- Update branch names given in compilation instructions
  (`doc/sphinx/develop/compiling/source-code.md`)

- Update entries in `doc/sphinx/_static/doc-versions.json`

- Update DOI in `README.rst`

- Make sure copyright years are correct:
  - `doc/sphinx/conf.py`
  - `License.txt`

- Add release notes to `doc/sphinx/reference/releasenotes/vX.Y.md`, linked from
  `doc/sphinx/reference/releasenotes/index.md`.
  - The [`graphql.py`](https://github.com/Cantera/cantera-release-guide/blob/main/graphql.py)
  and [`commiting-users.py`](https://github.com/Cantera/cantera-release-guide/blob/main/commiting-users.py)
  scripts in the `cantera-release-guide` repo can help with generating some of this
  content.

## After Merging Release Pull Request
- Once the CI completes, download the `docs` artifact, extract it, and add the contents
  to the `api-docs` repository.

- Update the recommended citation in `community.md` in the `cantera-website` repository,
  including the year and DOI.

- Create a tag for the revision corresponding to the release. If this is a feature
  release, also create a maintenance branch. Assuming `upstream` corresponds to the
  `Cantera/cantera` repository:
  ```
  git tag -a vX.Y.Z
  # Tag message: "Cantera X.Y.Z release"
  git push upstream vX.Y.Z
  git checkout -b X.Y
  git push upstream X.Y
  ```

- Create the release tarball. This requires a couple of steps to merge the
  `example_data` submodule in with the main source directories:
  ```sh
  # (after checking out the release commit)
  git archive --format=tar --prefix=cantera/ HEAD >source.tar
  cd data/example_data
  git archive --format=tar --prefix=cantera/data/example_data/ HEAD >../../example-data.tar
  cd ../..
  bsdtar czvf cantera-X.Y.Z.tar.gz @source.tar @example-data.tar
  rm source.tar example-data.tar
  ```
  - The trailing slash on the `--prefix` argument is _really important_. Otherwise,
  the files won't go into a folder.
  - Inspect the contents of the archive to make sure the directory structure looks correct:
    ```
    tar tvf cantera-X.Y.Z.tar.gz
    ```
  - The `bsdtar` tool can be installed on Debian/Ubuntu from the `libarchive-tools`
    package, or the `libarchive` Conda package. It is installed by default on macOS (and
    is the same as `tar`).

- Build officially-maintained packages
  - [PyPI](https://cantera.org/dev/develop/distribution-packages/pypi-sdist-wheel)
  - [conda-forge](https://cantera.org/dev/develop/distribution-packages/conda)
  - [Ubuntu](https://github.com/Cantera/cantera-ubuntu)

- Generate the detailed Changelog using the
  [`graphql.py`](https://github.com/Cantera/cantera-release-guide/blob/main/graphql.py)
  script from the `cantera-release-guide` repo and put the resulting list in the
  [wiki](https://github.com/Cantera/cantera/wiki)
