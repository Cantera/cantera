function a = isInlet(d)
% ISINLET  Determine whether a domain is an inlet.
% a = isInlet(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     1 if the domain is an inlet, and 0 otherwise.
%

a = int8(strcmp(domainType(d), 'inlet'));
end
