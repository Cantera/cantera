function h = enthalpy_mass(r)
% ENTHALPY_MASS - the specific enthalpy [J/kg].
%
%       h = enthalpy_mass(r)
%
%    returns the specific enthalpy of the reactor contents at the
%    end of the last call to 'advance' or 'step.'
%
%    See also: Reactor/intEnergy_mass, Reactor/entropy_mass,
%    Reactor/enthalpy_mole
%
h = reactormethods(27, reactor_hndl(r));
