function n = nElements(self)
% NELEMENTS  Get the number of elements in a mixture.
% n = nElements(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Number of elements in the input
%

n = mixturemethods(21, mix_hndl(self));
