function n = nSpecies(self)
% NSPECIES  Get the number of species in a mixture.
% n = nSpecies(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Number of species in the input
%

n = mixturemethods(24, mix_hndl(self));
