**************************
Frequently Asked Questions
**************************

Installation & Compilation
--------------------------

**How do I install Cantera on Windows?**

    Download the MSI installer for Cantera and the corresponding Python module
    from `SourceForge <https://sourceforge.net/projects/cantera/files/cantera/>`_.
    Choose between x86 and x64 based on the versions of Python and/or Matlab
    you want to work with. See :ref:`Windows Installation <sec-install-win>`
    for details.

**How do I install Cantera on Linux?**

    For Ubuntu, packages for the current stable version of Cantera are available
    in a PPA. See :ref:`Ubuntu Installation <sec-install-ubuntu>` for details.

    For other Linux distributions, download the source code (e.g.
    ``cantera-2.1.1.tar.gz``) from `SourceForge
    <https://sourceforge.net/projects/cantera/files/cantera/>`_ and follow the
    instructions in the :ref:`sec-compiling`.

**How do I install Cantera on Mac OS X?**

    Cantera can be installed using Homebrew. See :ref:`Mac OS X Installation
    <sec-install-osx>` for details.

**What do I do if compiling Cantera fails?**

    - Examine the output of the ``scons build`` command, especially anything
      identified as a       ``WARNING`` or ``ERROR``. Check for discrepancies
      with your expected configuration (e.g. not finding SUNDIALS even though
      you have it installed).
    - Check the contents of ``cantera.conf`` to make sure they are correct.
    - If any of the configuration tests (``Checking for...``) fail unexpectedly,
      look at the contents of ``config.log`` to determine the reason.
    - If none of these help identify the cause of the failure, consider asking
      for help on the Cantera Users' Group. If you decide to make a post, please
      include the following information:

      * The contents of ``cantera.conf`` and ``config.log``
      * The output of the ``scons build`` and ``scons build dump`` commands
        (you can direct this output to a file by running ``scons build >buildlog.txt 2>&1``)
      * The exact version of Cantera you are trying to compile, and how it was
        obtained (i.e. downloaded source tarball or the specific Git commit)
      * Your operating system, compiler versions, and the versions of any other
        relevant software.

**How do I debug issues with the SCons build system?**

    Sometimes, it is helpful to see all of the internal variables defined by
    SCons, either automatically or by the Cantera build scripts. To do this, add
    ``dump`` to your SCons command line. For example::

        $ scons build dump

    will show the variables that would be set during the ``build`` step. Note
    that in this case, the ``build`` step will not be executed.

    Alternatively, it is also possible to run SCons through the Python debugger, and set a breakpoint in the ``SConstruct`` file. For example::

        $ scons --debug=pdb build
        (Pdb) b /full/path/to/SConstruct:33
        (Pdb) cont

General
-------

**Which Cantera interface should I use?**

    If you're new to Cantera, the best interface to get started with is
    probably the Python interface. It offers most of the features of the
    C++ core in a much more flexible environment. Since all of the
    calculations are still done in C++, there is very little performance
    penalty to using the high-level language interfaces.

**Where can I find examples of how to use Cantera?**

    Cantera is distributed with many examples for the Python and Matlab
    interfaces, and a smaller number of examples for the C++ and Fortran
    interfaces. The Matlab, C++, and Fortran examples should be
    installed in the ``samples`` subdirectory of the Cantera installation
    directory, or they can be found in the ``samples`` subdirectory of the
    Cantera source directory.

    Examples for the Python interface can be found in the ``examples``
    subdirectory of the Cantera Python module installation directory, or in
    the ``interfaces/cython/cantera/examples`` subdirectory of the Cantera
    source directory.

**How should I cite Cantera?**

    The recommended citation for Cantera is as follows:

    David G. Goodwin, Harry K. Moffat, and Raymond L. Speth. *Cantera: An object-
    oriented software toolkit for chemical kinetics, thermodynamics, and
    transport processes*. http://www.cantera.org, 2017. Version 2.3.0.
    doi:10.5281/zenodo.170284

    The following BibTeX entry may also be used::

        @Misc{Cantera,
           author = "David G. Goodwin and Harry K. Moffat and Raymond L. Speth",
           title = "Cantera: An Object-oriented Software Toolkit for Chemical
                    Kinetics, Thermodynamics, and Transport Processes",
           year = 2017,
           note = "Version 2.3.0",
           howpublished = "\url{http://www.cantera.org}"
           doi = {10.5281/zenodo.170284}
        }

    If you are using a different version of Cantera, update the ``version`` and
    ``year`` fields accordingly.


Support and Bug Reporting
-------------------------

**What should I do if I think I've found a bug in Cantera?**

    - Check to see if you're using the most recent version of Cantera, and
      upgrade if not.
    - Check the `Issue Tracker
      <https://github.com/Cantera/cantera/issues>`_ to see if the issue
      has already been reported.
    - Try to generate a `minimal, complete, and verifiable example
      <http://stackoverflow.com/help/mcve>`_ that demonstrates the observed bug.
    - Create a new issue on the tracker. Include as much information as
      possible about your system configuration (operating system, compiler
      versions, Python versions, installation method, etc.)

**What information should I include in my bug report?**

    - The version of Cantera are you using, and how you installed it
    - The operating system you are using
    - If you compiled Cantera, what compiler you used, and what compilation
      options you specified
    - The version of Python or Matlab are you using, if applicable
    - The necessary *input* to generate the reported behavior
    - The full text of any error message you receive

**What should I do if I need help using Cantera?**

    You can join the `Cantera Users' Group
    <https://groups.google.com/forum/#!forum /cantera-users>`_ on Google
    Groups and ask a question there. Please use the search feature before
    posting to see if your question has been answered before. This group is
    moderated, so it may take some time for your posts to appear if you are a
    new member.
