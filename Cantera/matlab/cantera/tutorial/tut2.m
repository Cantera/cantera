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

gas = IdealGasMix('gri30.xml')

% Function 'IdealGasMix' constructs an object representing an ideal
% gas mixture by reading in attributes of the mixture from a file,
% which in this case is 'gri30.xml'. This file contains a complete
% specification of the GRI-Mech 3.0 reaction mechanism, including
% element data (name, atomic weight), species data (name, elemental
% composition, coefficients to compute thermodynamic and transport
% properties), and reaction data (stoichiometry, rate coefficient
% parameters). The file is written in an XML dialect understood by
% Cantera (CTML, or "Cantera Markup Language").
%
% Several reaction mechanism files in CTML format are included in the
% Cantera distribution, including ones that model high-temperature air
% and a hydrogen/oxygen reaction mechanism.  Under Windows, the
% installation program puts these files in 'C:\Program
% File\Common Files\Cantera.'  On a unix/linux machine, they are
% kept in the 'data/inputs' subdirectory within the Cantera
% installation directory.
%
%
% CK-format files
% ---------------
%
% Cantera also can use reaction mechanism files written in the format
% used in the Chemkin-II software package [1], which we will refer to
% as 'CK format'. Many gas-phase reaction mechanisms are available in
% this format. (See, for example, 
%  http://www.galcit.caltech.edu/EDL/mechanisms/library/library.html)

% To use a CK-format reaction mechanism file, call IdealGasMix with
% the file name, followed by the name of the thermodynamic property
% database (optional), followed by the name of the transport
% database (optional):
%
% gas1 = IdealGasMix('mech1.inp')
% gas2 = IdealGasMix('mech2.inp','therm.dat')
% gas3 = IdealGasMix('mech2.inp','therm.dat','tran.dat')
%
% If no transport database is specified, then transport property
% evaluation is disabled. 
%
% When function IdealGasMix is called with a CK file as an
% argument, an equivalent file in CTML format is written, and then
% the mixture attributes are read in from the CTML file. The next
% time you use the mechanism, you can simply specify the CTML file:
%
% gas2 = IdealGasMix('mech2.xml')
%
% Note that CTML files contain all required thermodynamic and
% transport data, and so no external databases need be specified.


% The Search Path
% ---------------
%
% Cantera looks in several standard directories for input files. As
% long as you have set CANTERA_ROOT correctly, you can use file
% 'gri30.xml' and the others that come with Cantera from any
% directory. But note that Cantera always looks in the local directory
% first. So if you have a modified file of the same name in the local
% directory, it will be used instead.
%
% If you have a directory where you keep files you want to use with
% Cantera, you can add it to the Cantera search path. For example, to
% add directory '/usr/local/lib/inputs' to the search path, call the
% MATLAB function 'addDirectory':
%
% addDirectory('/usr/local/lib/inputs')


%----------------------------------------------------------------
% [1] R. J. Kee, F. M. Rupley, and J. A. Miller, Sandia National
% Laboratories Report SAND89-8009 (1989).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end of tutorial 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

