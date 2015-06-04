function u = intEnergy_mass(r)
% INTENERGY_MASS  Get the specific internal energy.
% u = intEnergy_mass(r)
% See also: :mat:func:`enthalpy_mass`
%
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The specific internal energy of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     The internal energy is retrieved from the solution vector.
%     Units: J/kg
%

u = reactormethods(28, reactor_hndl(r));
