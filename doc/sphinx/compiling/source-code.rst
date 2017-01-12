
.. _sec-source-code:

Downloading the Cantera source code
===================================

Stable Release
--------------

* **Option 1**: Check out the code using Git::

    git clone --recursive https://github.com/Cantera/cantera.git
    cd cantera

  Then, check out the tag of the most recent stable version::

    git checkout tags/v2.3.0

  A list of all the tags can be shown by::

    git tag --list

* **Option 2**: Download the most recent source tarball from `Github
  <https://github.com/Cantera/cantera/releases>`_ and extract the
  contents.

Beta Release
------------

* Check out the code using Git::

    git clone --recursive https://github.com/Cantera/cantera.git
    cd cantera

  Then pick either **Option 1** or **Option 2** below.

* **Option 1**: Check out the tag with the most recent beta release::

    git checkout tags/v2.3.0b1

  Note that the most recent beta version might be older than the most recent
  stable release. A list of all the tags, including stable and beta versions can
  be shown by::

    git tag --list

* **Option 2**: Check out the branch with all the bug fixes leading to the
  next minor release of the stable version::

    git checkout 2.3

  This branch has all the work on the 2.3.x version of the software.

Development Version
-------------------

Check out the code using Git::

  git clone --recursive https://github.com/Cantera/cantera.git
  cd cantera

Note that by default, the ``master`` branch is checked out, containing all of
the feature updates and bug fixes to the code since the previous stable release.
The master branch is usually an "alpha" release, corresponding to the ``a`` in
the version number, and does not usually get a tag.
