function moles = phaseMoles(self, n)
% PHASEMOLES - moles of phase number 'n' (kmol).
%
if nargin == 2
    moles = mixturemethods(28, mix_hndl(self), n);
elseif nargin == 1
    np = nPhases(self);
    m = zeros(1,np);
    for n = 1:np
        m(n) = mixturemethods(28, mix_hndl(self), n);
    end
    moles = m;
else
    error('wrong number of arguments');
end
