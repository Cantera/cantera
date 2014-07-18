function moles = phaseMoles(self, n)
% PHASEMOLES  Get the number of moles of a phase in a mixture.
% moles = phaseMoles(self, n)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :param n:
%     Integer phase number in the input
% :return:
%     Moles of phase number ``n``. Units: kmol
%

if nargin == 2
    moles = mixturemethods(28, mix_hndl(self), n);
elseif nargin == 1
    np = nPhases(self);
    m = zeros(1, np);
    for n = 1:np
        m(n) = mixturemethods(28, mix_hndl(self), n);
    end
    moles = m;
else
    error('wrong number of arguments');
end
