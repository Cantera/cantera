function kr = revRateConstants(a)
%REVRATECONSTANTS Reverse reaction rate constants for all the reactions.
%
%    kr = revRateConstants(a)
%
%        Returns a column vector of the reverse rate constants of
%        all of the reactions.
%
kr = kinetics_get(a.id,16,0);
