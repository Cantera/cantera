function h = enthalpy_mass(r)
% ENTHALPY_MASS  The specific enthalpy of the reactor.
% h = enthalpy_mass(r)
%
% See also: :mat:func:`intEnergy_mass`
%
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The specific enthalpy of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     The enthalpy is retrieved from the solution vector. Units: J/kg
%

h = reactormethods(27, reactor_hndl(r));
