function setTemperature(self, T)
% SETTEMPERATURE  Set the mixture temperature.
% setTemperature(self, T)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :param T:
%     Temperature. Units: K
%

mixturemethods(5, mix_hndl(self), T);
