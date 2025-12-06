# Contributing to Cantera

As an open source project, Cantera relies on contributions from its users. We welcome
your participation in the project in any capacity that you are able to contribute,
whether that be through writing code, updating documentation, submitting bug reports, or
answering questions on the Users' Group. The following guidelines are intended to
provide a starting point for contributing code to Cantera.

* If you have an idea for significant change you'd like to implement, please create an
  issue on the
  [Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose)
  repository on GitHub or start a discussion on the
  [Cantera Users' Group](https://groups.google.com/g/cantera-users) to plan your
  modifications so that they can be implemented efficiently and in a way that doesn't
  conflict with any other planned future development
* If you are interested in contributing but don't have a particular goal in mind, there
  are a few ways to get started:
  * Look over the [Enhancements](https://github.com/Cantera/enhancements/issues) that are listed as "feature requests"
  * Look over the list of open [Issues](https://github.com/Cantera/cantera/issues)
  * Reach out to the [core developers](mailto:developers@cantera.org) to discuss areas
    where you could contribute effectively
* When making changes, write code that is consistent with the surrounding code
  (see the [Cantera style guidelines](https://cantera.org/dev/develop/style-guidelines.html))
* Add tests for any new features that you are implementing to either the
  GoogleTest-based test suite or the Python test suite. See
  [Writing Tests](https://cantera.org/dev/develop/writing-tests.html) for more
  information.
* Add examples that highlight new capabilities, or update existing examples to make use
  of new features.
* Cantera is licensed under a [BSD
  license](https://github.com/Cantera/cantera/blob/main/License.txt) which
  allows others to freely modify the code, and if your changes are accepted,
  then that code will be release under this license as well. The copyright for
  Cantera is held collectively by the contributors.
* You can find additional information about how Cantera is structured and tips for
  developing and debugging Cantera in the [Develop](https://cantera.org/dev/develop/)
  section of the Cantera website.
* Keep the eventual code review process in mind while you work. This article on
  [The (written) unwritten guide to pull requests](https://www.atlassian.com/blog/git/written-unwritten-guide-pull-requests) provides some good suggestions on how to
  prepare your contribution in a way that makes it efficient to review.

## Getting Credit for your Contributions

* Configure Git with your name and e-mail address before making any commits
  * From a terminal, run:
    ```shell
    git config --global user.name "Your Name"
    git config --global user.email "you@somewhere.org"
    ```
  * Make sure the e-mail you specify is linked to your GitHub account. You can view
    these settings on GitHub by clicking on your profile image in the upper right corner
    and selecting "Settings". Then, select the "Emails" section from the list on the
    left of the page.
* If you have made a significant contribution (more than just fixing a typo), please add
  your name, GitHub handle, and institutional affiliation (if desired) to the
  [Cantera authors](https://github.com/Cantera/cantera/blob/main/AUTHORS.md) file
  `AUTHORS.md` as part of your pull request.
* Do not add any author acknowledgements to individual files.

## Git Workflow

* Fork the `Cantera/cantera` repository on Github
* Clone your new repository or add it as a remote to an existing repository
* Check out the existing `main` branch, then start a new feature branch for your work.
  * Starting from the latest commit on the `main` branch, which is always the current
    development version, will minimize the likelihood of conflicts when merging your
    pull request.
  * Do not make commits directly to your `main` branch, or open a pull request from your
    `main` branch.
* As you make changes, commit them to your feature branch
  * Use descriptive commit messages (summary line of no more than 72 characters,
    followed by a blank line and a more detailed summary, if any)
  * Make related changes in a single commit, and unrelated changes in separate commits
  * Make sure that your commits do not include any undesired files, such as files
    produced as part of the build process or other temporary files.
  * Organize your commits and squash small "fixup" commits and reversions using Git's
    history-rewriting features such as `git rebase -i`; see
    <https://help.github.com/articles/about-git-rebase/>
  * *Do not merge your branch with `main`*. If your branch gets behind Cantera's `main`
    branch, you should *rebase* your branch onto the most recent `HEAD` commit of
    `main`. This allows us to keep a simple, linear project history.
  * Periodically run the test suite (`scons test`) to make sure that your
    changes are not causing any test failures.
* Push the changes on your new feature branch to your forked copy of the
  `Cantera/cantera` repository on GitHub.
* Submit a Pull Request on Github, from your forked copy. Check the results
  of the continuous-integration tests run using GitHub Actions and resolve
  any issues that arise.
* When responding to review comments on your pull request, please do not use the "accept
  suggestion" option, as this creates a very noisy version history.

The following pages describe the Git development workflows for some other projects. You
may also find some helpful advice from these projects:
* [Matplotlib](https://matplotlib.org/devdocs/devel/development_workflow.html)
* [SciPy](https://docs.scipy.org/doc/scipy/dev/contributor/development_workflow.html)

## Restrictions on Generative AI Usage

We expect authentic engagement in our community.

- Do not post output from Large Language Models or similar generative AI as comments on
  GitHub or our User Group, as such comments tend to be formulaic and low content.
- If you use generative AI tools as an aid in developing code or documentation changes,
  you are required to disclose the use of AI and ensure that you fully understand the
  proposed changes and can explain why they are the correct approach.

Make sure you have added value based on your personal competency to your contributions.
Just taking some input, feeding it to an AI and posting the result is not of value to
the project. To preserve precious core developer capacity, we reserve the right to
rigorously reject seemingly AI generated low-value contributions.
(_adopted from
[Matplotlib's AI Usage](https://matplotlib.org/devdocs/devel/contribute.html#restrictions-on-generative-ai-usage)
guidelines_)

## AI Disclosures and Authorship

Cantera is distributed under a BSD-3-Clause license with contributions from
California Institute of Technology, Sandia Corporation (with U.S. Government rights),
and the Cantera Developers. Because this license requires clear attribution and
provenance, contributors must follow the rules below when using generative AI tools.

### **1. AI Use Disclosure**

All pull requests must include a brief **AI Use Statement** in the PR description:

- **No generative AI was used.** This contribution was written entirely without AI
  assistance.
- **Limited use of generative AI.**
  Standard or boilerplate code snippets were generated with AI and manually reviewed;
  all design, logic, and implementation decisions were made by the contributor.
  Examples: IDE code-completions or brief LLM queries for common patterns.
- **Extensive use of generative AI.**
  Significant portions of code or documentation were generated with AI, including
  logic and implementation decisions. All generated code and documentation were
  reviewed and understood by the contributor. Examples: Output from agentic coding
  tools and/or substantial refactoring by LLMs (web-based or local).

**How to choose:**
- **Limited:** AI helped with syntax, boilerplate, or standard patterns; you designed
  the solution.
- **Extensive:** AI suggested the approach, algorithm, or structure; you validated it
  works.

### **2. Authorship and Licensing Requirements**

Contributors must ensure that:

- They have the rights to license their contribution under the project’s BSD-3-Clause
  terms.
- They do **not knowingly** submit AI-generated content derived from copyrighted or
  license-incompatible sources (examples: GPL code, proprietary libraries).
- All required attribution and copyright notices for contributed material are preserved
  as required under BSD-3-Clause.

### **3. Provenance Awareness**

As generative AI systems may emit material similar to third-party codebases,
contributors must exercise caution. Contributors should exercise reasonable care in:

- verifying that AI suggestions are sufficiently original,
- avoiding verbatim snippets from incompatible licenses, and
- ensuring that any third-party content incorporated intentionally follows proper
  attribution and is license-compatible.

### **4. Maintainer Review and Discretion**

Maintainers may:

- request clarification about the origin of AI-assisted code,
- ask the contributor to manually rewrite suspect portions, or
- reject contributions where provenance is unclear or license compatibility cannot be
  reasonably established.

These measures protect the integrity of the project’s licensing and long-term
maintainability.
