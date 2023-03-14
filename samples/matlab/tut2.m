% Tutorial 2: Working with input files
%
%   Topics:
%     - using functions 'Solution' and 'importInterface'
%     - input files distributed with Cantera
%     - the Cantera search path
%     - CTML files
%     - converting from CK format
%
% Keywords: tutorial, input files

help tut2

t0 = cputime;

% In the last tutorial, we used function GRI30 to create an object
% that models an ideal gas mixture with the species and reactions of
% GRI-Mech 3.0. Another way to do this is shown here, with statements
% added to measure how long this takes:

gas1 = Solution('gri30.yaml', 'gas', 'mixture-averaged');
msg = sprintf('time to create gas1: %f', cputime - t0)

% Function 'Solution' constructs an object representing a phase of
% matter by reading in attributes of the phase from a file, which in
% this case is 'gri30.yaml'. This file contains a complete specification
% of the GRI-Mech 3.0 reaction mechanism, including element data
% (name, atomic weight), species data (name, elemental composition,
% coefficients to compute thermodynamic and transport properties), and
% reaction data (stoichiometry, rate coefficient parameters). The file
% is written in a format understood by Cantera, which is described in
% the document "Defining Phases and Interfaces."

% On some systems, processing long files can be a little slow. Loading
% a file again takes much less time, because Cantera 'remembers' files
% it has already processed and doesn't need to read them in again:

t0 = cputime;
gas1b = Solution('gri30.yaml', 'gas', 'none');
msg = sprintf('time to create gas1b: %f', cputime - t0)

% Since GRI-Mech is a rather small mechanism, there might not be much
% difference in these times.

% YAML files distributed with Cantera
%------------------------------------

% Several reaction mechanism files in this format are included in the
% Cantera distribution, including ones that model high-temperature
% air, a hydrogen/oxygen reaction mechanism, and a few surface
% reaction mechanisms. Under Windows, these files may be located in
% 'C:\Program Files\Common Files\Cantera', or in 'C:\cantera\data',
% depending on how you installed Cantera and the options you
% specified. On a Unix/Linux/macOS machine, they are usually kept
% in the 'data' subdirectory within the Cantera installation
% directory.

% If for some reason Cantera has difficulty finding where these files
% are on your system, set environment variable CANTERA_DATA to the
% directory where they are located. Alternatively, you can call function
% adddir to add a directory to the Cantera search path:
adddir('/usr/local/cantera/my_data_files');

% Cantera input files are plain text files, and can be created with
% any text editor. See the document 'Defining Phases and Interfaces'
% for more information.


% Importing multiple phases or interfaces
% ---------------------------------------

% A Cantera input file may contain more than one phase specification,
% or may contain specifications of interfaces (surfaces). Here we
% import definitions of two bulk phases and the interface between them
% from file diamond.yaml:

gas2 = Solution('diamond.yaml', 'gas');        % a gas

diamond = Solution('diamond.yaml','diamond');  % bulk diamond

diamonnd_surf = importInterface('diamond.yaml','diamond_100',...
                                gas2, diamond);

% Note that the bulk (that is, 3D) phases that participate in the surface
% reactions must also be passed as arguments to importInterface.

% Converting CK-format files
% --------------------------

% Many existing reaction mechanism files are in "CK format," by which
% we mean the input file format developed for use with the Chemkin-II
% software package. [See R. J. Kee, F. M. Rupley, and J. A. Miller,
% Sandia National Laboratories Report SAND89-8009 (1989).]

% Cantera comes with a converter utility program 'ck2yaml' (or
% 'ck2yaml.exe') that converts CK format into YAML format. This
% program should be run from the command line first to convert any CK
% files you plan to use into YAML format.
%
% Here's an example of how to use it:
%
% ck2yaml --input=mech.inp --thermo=therm.dat --transport=tran.dat --name=mymech

clear all
cleanup
