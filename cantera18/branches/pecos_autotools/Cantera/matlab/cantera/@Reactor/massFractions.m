function y = massFractions(r)
% MASSFRACTION - Mass fractions of reactor contents after last call
% to 'advance' or 'step'. 
%   
nsp = nSpecies(r.contents);
ir = reactor_hndl(r);
for k = 1:nsp
  yy(k) = reactormethods(30, ir, k-1);
end
y = yy;

