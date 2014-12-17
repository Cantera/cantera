function setThermalResistance(w, r)
% SETTHERMALRESISTANCE  Set the thermal resistance.
% setThermalResistance(w, r)
% :param w:
%     Instance of class :mat:func:`Wall`
% :param r:
%     Double, thermal resistance. Units: K*m**2/W
%

wallmethods(6, wall_hndl(w), r)
