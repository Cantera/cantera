function mu = chemPotentials(self)
% CHEMPOTENTIALS  Get the chemical potentials of species in a mixture.
% mu = chemPotentials(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :return:
%     Vector of chemical potentials. Units: J/kmol
%

mu = mixturemethods(41, mix_hndl(self));
