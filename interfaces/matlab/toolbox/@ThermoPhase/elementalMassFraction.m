function elMassFrac = elementalMassFraction(tp, element)
% ELEMENTALMASSFRACTION  Determine the elemental mass fraction in gas object.
% elMassFrac = elementalMassFraction(tp, element)
% :param tp:
%     Object representing the gas, instance of class :mat:func:`Solution`,
%     and an ideal gas. The state of this object should be set to an
%     estimate of the gas state before calling elementalMassFraction.
% :param element:
%     String representing the element name.
% :return:
%     Elemental mass fraction within a gas object.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input parameters
%

if nargin ~= 2
    error('elementalMassFraction expects two input arguments.');
end
if ~isIdealGas(tp)
    error('Gas object must represent an ideal gas mixture.');
end
if ~ischar(element)
    error('Element name must be of format character.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the elemental mass fraction in a gas object. Equation used is
% elMassFrac = sum of nAtoms(k,m)*Mel(m)*Y(k)/mw(k) where nAtoms(k,m) is
% the number of atoms of element, m, in species, k; Mel(m) is the atomic
% weight of the element, m; Y(k) is the mass fraction of species,k, in the
% gas object; and mw(k) is the molecular weight of species, k.
%

n = nSpecies(tp);
massFrac = massFractions(tp);
spec = speciesNames(tp);
eli = elementIndex(tp, element);
M = atomicMasses(tp);
Mel = M(eli);
MW = molecularWeights(tp);
% Initialize the element mass fraction as zero.
elMassFrac = 0.0;
% Use loop to perform summation of elemental mass fraction over all species.
for i = 1:n
    natoms(i) = nAtoms(tp,spec{i},element);
    mw(i) = MW(i);
    Y(i) = massFraction(tp,spec{i});
    elMassFrac = elMassFrac + (natoms(i)*Mel*Y(i))/mw(i);
end
end
