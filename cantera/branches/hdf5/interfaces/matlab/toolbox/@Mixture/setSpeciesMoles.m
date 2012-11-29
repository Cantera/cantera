function setSpeciesMoles(self, moles)
% SETSPECIESMOLES - Set the moles of the species [kmol]. The moles may
% be specified either as a string, or as an array. If an array is
% used, it must be dimensioned at least as large as the total number
% of species in the mixture. Note that the species may belong to any
% phase, and unspecified species are set to zero.
%
%     >> setSpeciesMoles(mix, 'C(s):1.0, CH4:2.0, O2:0.2');
%
mixturemethods(8, mix_hndl(self), moles);
