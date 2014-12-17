function setHeatTransferCoeff(w, u)
% SETHEATTRANSFERCOEFF  Set the heat transfer coefficient.
% setHeatTransferCoeff(w, u)
% :param w:
%     Instance of class :mat:func:`Wall`
% :param u:
%     Heat transfer coefficient of the wall. Units: W/(m**2-K)
%

wallmethods(7, wall_hndl(w), u);
