# Distributing Conda Package via conda-forge

The Conda recipe hosted at <https://github.com/conda-forge/cantera-feedstock> is used
to build packages for both stable and pre-release versions of Cantera.

## Stable Cantera releases

- fork the repository above
- Update `recipe/meta.yaml`:
  - Change the header line that looks like `{% set version = "3.1.0" %}`
  - In the `source` section, update the `url` and `sha256` for both the main `cantera`
    repo as well as the `cantera-example-data` repo.
    - The hash in the `url` is the hash of the commit corresponding to the tagged
      release.
    - The easiest way to get the `sha256` of the archive is to download it from that URL
      and run `sha256sum` on that file.
  - Update the value of `CT_GIT_COMMIT` with the same hash as used in the URL for the
    `cantera` repository.
  - Under `build`, reset the `number` to 0.
- Do whatever else is needed to bring the conda-forge recipe up to date, for example
  by using conda-smithy.
- Create a pull request against the `main` branch of `conda-forge/cantera-feedstock`
- Once the PR is merged, you should be able to see the status of the final package
  builds in [Azure DevOps](https://dev.azure.com/conda-forge/feedstock-builds/_build?definitionId=11466&_a=summary)

## Cantera pre-releases

- Follow the same steps as above but create the pull request against the `dev` branch of
  `conda-forge/cantera-feedstock`
- Whether your starting point should be the existing `main` or `dev` branch is a bit
  tricky. In either case, there are likely to be conflicts to work out. Use caution to
  make sure that important updates to the build aren't lost when merging branches.
