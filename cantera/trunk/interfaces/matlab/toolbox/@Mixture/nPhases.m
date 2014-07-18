function n = nPhases(self)
% NPHASES  Get the number of phases in a mixture.
% n = nPhases(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Number of phases in the input
%

n = mixturemethods(19, mix_hndl(self));
