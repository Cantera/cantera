function y = massFractions(r)
% MASSFRACTIONS  Get the mass fractions of the species.
% y = massFractions(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The mass fractions of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%

nsp = nSpecies(r.contents);
ir = reactor_hndl(r);
y = zeros(1, nsp);
for k = 1:nsp
    y(k) = reactormethods(30, ir, k-1);
end
