function setPressure(self, P)
% SETPRESSURE  Set the pressure of the mixture.
% setPressure(self, P)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :param P:
%     Pressure. Units: Pa
%

mixturemethods(6, mix_hndl(self), P);
