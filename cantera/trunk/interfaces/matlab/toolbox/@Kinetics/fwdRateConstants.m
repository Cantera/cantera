function kf = fwdRateConstants(a)
%FWDRATECONSTANTS Forward reaction rate constants for all the reactions.
%
%    kf = fwdRateConstants(a)
%
%        Returns a column vector of the forward rate constants of
%        all of the reactions.
%
kf = kinetics_get(a.id,15,0);
