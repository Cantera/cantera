%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Tutorial 2: Using your own reaction mechanism files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Function 'IdealGasMix'
% ----------------------

% In the last tutorial, we used function GRI30 to create an object
% that models an ideal gas mixture with the species and reactions of
% GRI-Mech 3.0. Another way to do this is shown here:

gas = IdealGasMix('gri30.cti')

% Function 'IdealGasMix' constructs an object representing an ideal
% gas mixture by reading in attributes of the mixture from a file,
% which in this case is 'gri30.cti'. This file contains a complete
% specification of the GRI-Mech 3.0 reaction mechanism, including
% element data (name, atomic weight), species data (name, elemental
% composition, coefficients to compute thermodynamic and transport
% properties), and reaction data (stoichiometry, rate coefficient
% parameters). The file is written in a format understood by
% Cantera, which is described in the document "Defining Phases and
% Interfaces."
%
% Several reaction mechanism files in this format are included in the
% Cantera distribution, including ones that model high-temperature air
% and a hydrogen/oxygen reaction mechanism.  Under Windows, the
% installation program puts these files in 'C:\Program
% File\Common Files\Cantera.'  On a unix/linux machine, they are
% kept in the 'data' subdirectory within the Cantera
% installation directory.
%
%
% CK-format files
% ---------------
%
% Cantera also comes with a converter utility for reaction mechanism
% files written in the format used in the Chemkin-II software package
% [1], which we will refer to as 'CK format'. Many gas-phase reaction
% mechanisms are available in this format. (See, for example,
% http://www.galcit.caltech.edu/EDL/mechanisms/library/library.html)

% To use a CK-format reaction mechanism file, from the command line 
% type:
%
% ck2cti -i mech.ck -t therm.dat -tr tran.dat -id mymechname > mech.cti
%
% Here therm.dat is a CK-format file containing species thermo
% data, and tran.dat is a Chemkin-compatible transport database. If
% transport properties are not needed, the transport database can
% be omitted, and if all species thermo data are in the mechanism
% file, the thermo database can also be omitted.

% How does Cantera find .cti input files?  Cantera always looks in the
% local directory first. If it is not there, Cantera looks for it on
% its search path. It looks for it in the data directory specified
% when Cantera was built (by default this is /usr/local/cantera/data
% on unix systems). If you define environment variable
% CANTERA_DATA_DIR, it will also look there, or else you can call
% function addDirectory to add a directory to the search path.

% Warning: when Cantera reads a .cti input file, wherever it is
% located, it always writes a file of the same name but with extension
% .xml *in the local directory*. If you happen to have some other file
% by that name, it will be overwritten. Once the XML file is created,
% you can use it instead of the .cti file, which will result in
% somewhat faster startup.




%----------------------------------------------------------------
% [1] R. J. Kee, F. M. Rupley, and J. A. Miller, Sandia National
% Laboratories Report SAND89-8009 (1989).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end of tutorial 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

