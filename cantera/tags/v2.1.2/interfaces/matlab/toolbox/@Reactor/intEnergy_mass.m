function u = intEnergy_mass(r)
% INTENERGY_MASS - the specific internal energy [J/kg].
%
%       u = intEnergy_mass(r)
%
%    returns the specific internal energy of the reactor contents at
%    the end of the last call to 'advance' or 'step.'
%
%    See also: Reactor/enthalpy_mass, Reactor/entropy_mass,
%    Reactor/enthalpy_mole
%
u = reactormethods(28, reactor_hndl(r));
