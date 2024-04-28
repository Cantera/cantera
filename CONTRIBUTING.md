# Contributing to Cantera

* For significant changes, please create an issue on the
  [Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose)
  repository on GitHub or start a discussion on the Cantera Users' Group to plan your
  modifications so that they can be implemented efficiently and in a way that doesn't
  conflict with any other planned future development
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
  your name, GitHub handle, and institutional affiliation (if desired) to the `AUTHORS`
  file as part of your pull request.
* Do not add any author acknowledgements to individual files.

## Git Workflow

* Fork the `Cantera/cantera` repository on Github
* Clone your new repository or add it as a remote to an existing repository
* Check out the existing `main` branch, then start a new feature branch for your work
* As you make changes, commit them to your feature branch
  * Use descriptive commit messages (summary line of no more than 72 characters,
    followed by a blank line and a more detailed summary, if any)
  * Make related changes in a single commit, and unrelated changes in separate commits
  * Make sure that your commits do not include any undesired files, such as files
    produced as part of the build process or other temporary files.
  * Use Git's history-rewriting features (such as `git rebase -i`; see
    <https://help.github.com/articles/about-git-rebase/>) to organize your commits
    and squash "fixup" commits and reversions.
  * Do not merge your branch with `main`. If needed, you should rebase your branch
    onto the most recent `HEAD` commit of `main`.
  * Periodically run the test suite (`scons test`) to make sure that your
    changes are not causing any test failures.
* Push the changes on your new feature branch to your forked copy of the
  `Cantera/cantera` repository on GitHub.
* Submit a Pull Request on Github, from your forked copy. Check the results
  of the continuous-integration tests run using GitHub Actions and resolve
  any issues that arise.

Additional discussion of good Git & Github workflow is provided at
<http://matplotlib.org/devel/gitwash/development_workflow.html> and
<https://docs.scipy.org/doc/numpy-1.15.0/dev/gitwash/development_workflow.html>
