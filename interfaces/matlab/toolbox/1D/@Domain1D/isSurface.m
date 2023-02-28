function a = isSurface(d)
% ISSURFACE  Determine if a domain is a surface.
% a = isSurface(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     1 if the domain is a surface, and 0 otherwise.
%

a = int8(strcmp(domainType(d), 'surface'));
end
