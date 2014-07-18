function setProfile(d, n, p)
% SETPROFILE  Set the profile of a component.
% d = setProfile(d, n, p)
% Convenience function to allow an instance of :mat:func:`Domain1D` to
% have a profile of its components set when it is part of a :mat:func:`Stack`.
%
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param n:
%     Integer index of component, vector of component indices, string
%     of component name, or cell array of strings of component names.
% :param p:
%     n x 2 array, whose columns are the relative (normalized) positions
%     and the component values at those points. The number of positions
%     ``n`` is arbitrary.
%

if d.stack == 0
    error('Install domain in stack before calling setProfile.');
end

setProfile(d.stack,domainIndex(d), n, p);
