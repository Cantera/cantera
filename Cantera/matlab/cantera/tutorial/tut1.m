% Tutorial 1:  Getting started
%   
%  Topics:
%    - creating a gas mixture
%    - setting the state
%    - cleaning up
%
help tut1

% Start MATLAB, and at the prompt type:

gas1 = GRI30

% If you have successfully installed the Cantera toolbox, 
% you should see something like this:
%
%
%       temperature                  300  K
%       pressure                  101325  Pa
%       density                 0.081896  kg/m^3
%       mean mol. weight         2.01594  amu
%
%                           X                 Y   
%                     -------------     ------------ 
%                H2   1.000000e+000     1.000000e+000
%                 H   0.000000e+000     0.000000e+000
%                 O   0.000000e+000     0.000000e+000
%                O2   0.000000e+000     0.000000e+000
%                OH   0.000000e+000     0.000000e+000
%               H2O   0.000000e+000     0.000000e+000
%               HO2   0.000000e+000     0.000000e+000
%              H2O2   0.000000e+000     0.000000e+000
%                 C   0.000000e+000     0.000000e+000
%                CH   0.000000e+000     0.000000e+000
%               CH2   0.000000e+000     0.000000e+000
%            CH2(S)   0.000000e+000     0.000000e+000
%               CH3   0.000000e+000     0.000000e+000
%               CH4   0.000000e+000     0.000000e+000
%                CO   0.000000e+000     0.000000e+000
%               CO2   0.000000e+000     0.000000e+000
%               HCO   0.000000e+000     0.000000e+000
%              CH2O   0.000000e+000     0.000000e+000
%             CH2OH   0.000000e+000     0.000000e+000
%              CH3O   0.000000e+000     0.000000e+000
%             CH3OH   0.000000e+000     0.000000e+000
%               C2H   0.000000e+000     0.000000e+000
%              C2H2   0.000000e+000     0.000000e+000
%              C2H3   0.000000e+000     0.000000e+000
%              C2H4   0.000000e+000     0.000000e+000
%              C2H5   0.000000e+000     0.000000e+000
%              C2H6   0.000000e+000     0.000000e+000
%              HCCO   0.000000e+000     0.000000e+000
%             CH2CO   0.000000e+000     0.000000e+000
%             HCCOH   0.000000e+000     0.000000e+000
%                 N   0.000000e+000     0.000000e+000
%                NH   0.000000e+000     0.000000e+000
%               NH2   0.000000e+000     0.000000e+000
%               NH3   0.000000e+000     0.000000e+000
%               NNH   0.000000e+000     0.000000e+000
%                NO   0.000000e+000     0.000000e+000
%               NO2   0.000000e+000     0.000000e+000
%               N2O   0.000000e+000     0.000000e+000
%               HNO   0.000000e+000     0.000000e+000
%                CN   0.000000e+000     0.000000e+000
%               HCN   0.000000e+000     0.000000e+000
%              H2CN   0.000000e+000     0.000000e+000
%              HCNN   0.000000e+000     0.000000e+000
%              HCNO   0.000000e+000     0.000000e+000
%              HOCN   0.000000e+000     0.000000e+000
%              HNCO   0.000000e+000     0.000000e+000
%               NCO   0.000000e+000     0.000000e+000
%                N2   0.000000e+000     0.000000e+000
%                AR   0.000000e+000     0.000000e+000
%              C3H7   0.000000e+000     0.000000e+000
%              C3H8   0.000000e+000     0.000000e+000
%            CH2CHO   0.000000e+000     0.000000e+000
%            CH3CHO   0.000000e+000     0.000000e+000
%
% What you have just done is to create an object ("gas1") that
% implements GRI-Mech 3.0, the 53-species, 325-reaction natural gas
% combustion mechanism developed by Gregory P. Smith, David M. Golden,
% Michael Frenklach, Nigel W. Moriarty, Boris Eiteneer, Mikhail
% Goldenberg, C. Thomas Bowman, Ronald K. Hanson, Soonho Song, William
% C. Gardiner, Jr., Vitali V. Lissianski, and Zhiwei Qin. (See
% http://www.me.berkeley.edu/gri_mech/ for more information about
% GRI-Mech 3.0.)
%
% The object created by GI30 has properties you would expect for a gas
% mixture - it has a temperature, a pressure, species mole and mass
% fractions, etc. As we'll soon see, it has many other properties too.
%
% The summary of the state of 'gas1' printed above shows that new
% objects created by function GRI30 start out with a temperature of
% 300 K, a pressure of 1 atm, and have a composition that consists of
% only one species, in this case hydrogen. There is nothing special
% about H2 - it just happens to be the first species listed in the
% input file defining GRI-Mech 3.0 that the 'GRI30' function reads. In
% general, the species listed first will initially have a mole
% fraction of 1.0, and all of the others will be zero.

%  Setting the state
%  -----------------

% The state of the object can easily be changed. For example,

setTemperature(gas1, 1200)

% sets the temperature to 1200 K. (Cantera always uses SI units.)

% Notice in the summary of properties that MATLAB prints after this
% command is executed that the temperature has been changed as
% requested, but the pressure has changed too. The density and
% composition have not. 
%
% When setting properties individually, some convention needs to be
% adopted to specify which other properties are held constant. This is
% because thermodynamics requires that *two* properties (not one) in
% addition to composition information be specified to fix the
% intensive state of a substance (or mixture).
%
% Cantera adopts the following convention: only one of the set
% (temperature, density, mass fractions) is altered by setting any
% single property. In particular:
%
% a) Setting the temperature is done holding density and
%    composition  fixed. (The pressure changes.)

% b) Setting the pressure is done holding temperature and
%    composition fixed. (The density changes.)
% 
% c) Setting the composition is done holding temperature
%    and density fixed. (The pressure changes).
%


% Setting multiple properties: the 'set' method
% ---------------------------------------------

% If you want to set multiple properties at once, use the 'set'
% method. (Note: a 'method' is just the term for a function that acts
% on an object. In MATLAB, methods take the object as the first
% argument.)

set(gas1, 'Temperature', 900.0, 'Pressure', 1.e5)

% This statement sets both temperature and pressure at the same
% time. Any number of property/value pairs can be specified in a
% call to 'set'. For example, the following sets the mole fractions
% too:
set(gas1, 'Temperature', 900.0, 'Pressure', 1.e5, 'MoleFractions', ...
	  'CH4:1,O2:2,N2:7.52')

% The 'set' method also accepts abbreviated property names:

set(gas1,'T',900.0,'P',1.e5,'X','CH4:1,O2:2,N2:7.52')

% Either version results in
%
%       temperature                  900  K
%       pressure                  100000  Pa
%       density                   0.3693  kg/m^3
%       mean mol. weight         27.6332  amu
%
%                           X                 Y   
%                     -------------     ------------ 
%                O2   1.901141e-001     2.201489e-001
%               CH4   9.505703e-002     5.518732e-002
%                N2   7.148289e-001     7.246638e-001
%

% Other properties may also be set using 'set', including some that
% can't be set individually. The following property pairs may be
% set: (Enthalpy, Pressure), (IntEnergy, Volume), (Entropy,
% Volume), (Entropy, Pressure). In each case, the values of the
% extensive properties must be entered *per unit mass*. 

% Setting the enthalpy and pressure:
set(gas1, 'Enthalpy', 2*enthalpy_mass(gas1), 'Pressure', 2*oneatm);


% The composition above was specified using a string. The format is a
% comma-separated list of <species name>:<relative mole numbers>
% pairs. The mole numbers will be normalized to produce the mole
% fractions, and therefore they are 'relative' mole numbers.  Mass
% fractions can be set in this way too by changing 'X' to 'Y' in the
% above statement.

% The composition can also be set using an array, which can be
% either a column vector or a row vector but must have the same
% size as the number of species. For example, to set all 53 mole
% fractions to the same value, do this:

x = ones(53,1);   % a column vector of 53 ones
set(gas1, 'X', x)

% To set the mass fractions to equal values:
set(gas1, 'Y', x)


% This clears all Matlab objects created
clear all

% and this clears all Cantera objects created
cleanup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end of tutorial 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




