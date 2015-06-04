function t = temperature(self)
% TEMPERATURE  Get the temperature of a mixture.
% t = temperature(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Temperature (K)
%

t = mixturemethods(25, mix_hndl(self));
