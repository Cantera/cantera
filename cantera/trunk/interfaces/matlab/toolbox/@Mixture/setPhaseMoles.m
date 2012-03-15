function setPhaseMoles(self, n, moles)
% SETPHASEMOLES - set the number of moles of phase number 'n' to
% 'moles' (kmol).
%
mixturemethods(7, mix_hndl(self), n, moles);
