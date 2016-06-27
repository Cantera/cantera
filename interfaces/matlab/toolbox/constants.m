function [atm, r] = constants
% CONSTANTS  Get the values of important constants.
% [atm,r] = constants
% Deprecated. To be removed after Cantera 2.3. Use :mat:func:`oneatm` and
% :mat:func:`gasconstant` instead.
%
% :return:
%     If one output argument is given, returns one atmosphere in
%     Pascals. If two output arguments are given, returns one
%     atmosphere in Pascals and the universal gas constant in
%     J/kmol-K.
%

warning('This function is deprecated and will be removed after Cantera 2.3.');
atm = 101325.0;
r = 8314.4621;
