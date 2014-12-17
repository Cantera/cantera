function n = speciesIndex(self, k, p)
% SPECIESINDEX  Get the index of a species in a mixture.
% n = speciesIndex(self, k, p)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :param name:
%     Name of the speces whose index is desired
% :return:
%     Index of species with name ``name``
%

n = mixturemethods(23, mix_hndl(self), k, p);
