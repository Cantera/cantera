function mu = chemPotentials(self)
% CHEMPOTENTIALS - Chemical potentials of all species in all phases
mu = mixturemethods(41, mix_hndl(self));
