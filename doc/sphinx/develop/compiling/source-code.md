(sec-source-code)=
# Downloading the Cantera source code

## Stable Release

- **Option 1**: Check out the code using Git (recommended):

  ```bash
  git clone --recursive https://github.com/Cantera/cantera.git
  cd cantera
  ```

  Then, check out the tag of the most recent stable version:

  ```bash
  git checkout tags/v3.1.0
  git submodule update
  ```

  A list of all the tags can be shown by:

  ```bash
  git tag --list
  ```

- **Option 2**: Download the most recent source `.tar.gz` or `.zip` file from
  [GitHub](https://github.com/Cantera/cantera/releases) and extract the contents. In
  this case, several dependencies that are linked to the Cantera Git repository will not
  be available and will need to be installed elsewhere on your system.

```{button-ref} configure-build
:color: primary
:shadow:
:align: right
Next: Configure & Build Cantera
```

## Beta Release

- Check out the code using Git:

  ```bash
  git clone --recursive https://github.com/Cantera/cantera.git
  cd cantera
  ```

  Then pick either **Option 1** or **Option 2** below.

- **Option 1**: Check out the tag with the most recent beta release:

  ```bash
  git checkout tags/v3.1.0b1
  git submodule update
  ```

  Note that the most recent beta version might be older than the most recent stable
  release. A list of all the tags, including stable and beta versions can be shown by:

  ```bash
  git tag --list
  ```

- **Option 2**: Check out the branch with all the bug fixes leading to the next minor
  release of the stable version:

  ```bash
  git checkout 3.1
  git submodule update
  ```

  This branch has all the work on the 3.1.x version of the software.

  If you've already checked out the 3.1 branch, you can get the latest updates from the
  main Cantera repository and synchronize your local repository by running:

  ```bash
  git checkout 3.1
  git fetch --all
  git pull --ff-only
  ```

```{button-ref} configure-build
:color: primary
:shadow:
:align: right
Next: Configure & Build Cantera
```

## Development Version

Check out the code using Git:

```bash
git clone --recursive https://github.com/Cantera/cantera.git
cd cantera
```

Note that by default, the `main` branch is checked out, containing all of the feature
updates and bug fixes to the code since the previous stable release. The main branch is
usually an alpha release, corresponding to the `a` in the version number, and does not
usually get a tag.

If you've previously checked out the repository, and haven't made any changes locally,
you can get the latest updates from the main Cantera repository and synchronize your
local repository by running:

```bash
git checkout main
git fetch --all
git pull --ff-only
```

```{button-ref} configure-build
:color: primary
:shadow:
:align: right
Next: Configure & Build Cantera
```
