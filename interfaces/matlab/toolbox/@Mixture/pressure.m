function p = pressure(self)
% PRESSURE  Get the pressure of the mixture.
% p = pressure(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Pressure. Units: Pa
%

p = mixturemethods(26, mix_hndl(self));
