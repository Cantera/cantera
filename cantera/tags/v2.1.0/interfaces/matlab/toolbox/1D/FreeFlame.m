function m = FreeFlame(gas, id)
% FREEFLAME - Freely-propagating flat flame
%
%    Return a Domain1D instance representing a freely-propagating
%    adiabatic flame
%
m = Domain1D(1, gas, 2);
if nargin == 1
    setID(m,'flame');
else
    setID(m,id);
end
