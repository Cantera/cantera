# Continuous Integration

Cantera uses [GitHub Actions](https://docs.github.com/en/actions) to run a suite of tests and code quality checks to help avoid regressions and make sure Cantera's code has a consistent style. Many of these checks run automatically on each push to a pull request and a few run only when a pull request is merged into the `main` branch. This section describes how the workflows are set up and requirements for adding or changing checks.

## Existing Workflows

All the workflows that run on GitHub Actions are in the `.github/workflows` folder:

- `linters.yml`: Has actions that check that code is formatted correctly. Runs on every push for a pull request, and every commit to `main`.
- `main.yml`: Has actions that run the test suite to avoid unexpected changes in behavior. Runs on every push for a pull request, and every commit to `main`.
- `packaging.yml`: Has actions that check our packaging workflows. Runs on every commit to `main`.
- `post-merge-tests.yml`: Has actions that take a long time to run, or that run the test suite on less commonly-used platforms. Runs on every commit to `main`.
- `zizmor.yml`: Has an action to run [zizmor](https://woodruffw.github.io/zizmor/), a program that statically checks workflow files for best-practices. Runs on every push for a pull request, and every commit to `main`.

## Modifying or Adding Workflows

Workflows can be modified to include additional jobs or modify the behavior of existing jobs. Workflows can also be added if another category of check becomes useful. The best bet is to follow the existing patterns in the workflow files while making updates as this will help make sure that modified jobs work properly. Keep in mind the following tips as well:

1. Steps in a job that use a 3rd party action with the `uses` keyword must pin the action to a SHA hash. The line should be formatted as:

   ```yaml
   - name: Step
     uses: <organization>/<repository>@<40 character commit hash> # <version tag>
   ...
   ```

   Pinning actions in this way makes sure that our actions are always using the same version of a the 3rd party code that is expected.

   1. Steps developed by GitHub in the `actions` organization do not need to be pinned to a SHA hash but they can be if desired. If a SHA has is not used, steps from the `actions` repository _must_ be pinned to a tag or branch that matches the regex pattern `^v*$`.
   2. The versions of steps are kept up-to-date by Dependabot, an automated service that is configured to run monthly. Dependabot will update the hash and the version tag whenever updates for an action are available.
2. When checking out a repository in a job, please set `persist-credentials: false` in the `with` block. For example:

   ```yaml
   - name: Checkout the repository
     uses: actions/checkout@v4
     with:
       persist-credentials: false
       <other options>
    ```
3. Workflows and jobs should be configured with the [minimal permissions](https://docs.github.com/en/actions/writing-workflows/choosing-what-your-workflow-does/controlling-permissions-for-github_token) necessary to complete their steps. In general, this means you should add the following configuration at the top-level of every workflow:

   ```yaml
   permissions: {}
   ```

   If jobs need additional permissions, they should be configured per-job to make sure the extra permissions are not given to other jobs in the workflow as in the following example:

   ```yaml
   jobs:
     job-name:
       permissions:
         contents: read
         ...
   ```
4. Do not use template substitution with `${{ }}` directly in scripts. Doing so can result in unexpected code being executed if the template string contains something that Bash can interpret. Instead, assign the template to an environment variable and use the variable in the script:

   ```yaml
   # DO NOT DO THIS
   - run: echo "${{ github.ref }}"

   # DO THIS INSTEAD
   - run: echo "$GH_REF"
     env:
       GH_REF: "${{ github.ref }}"
   ```

   For more information, see the [Zizmor documentation](https://woodruffw.github.io/zizmor/audits/#template-injection) about template injection.
5. You can run Zizmor locally to statically check changes in workflow files with the following command:

   ```bash
   GH_TOKEN=$(gh auth token) zizmor .github/workflows
   ```

   This uses the [`gh` CLI program](https://cli.github.com/manual/) to get a token to run online checks with Zizmor.
