function setMoleFractions(d, x)
% SETMOLEFRACTIONS  Set the mole fractions.
% d = setMoleFractions(d, x)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param x:
%     String specifying the species and mole fractions in
%     the format ``'SPEC:X,SPEC2:X2'``.
%

domain_methods(d.dom_id, 62, x);
