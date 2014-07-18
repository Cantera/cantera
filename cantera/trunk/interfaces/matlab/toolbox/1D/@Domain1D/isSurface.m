function a = isSurface(d)
% ISSURFACE  Determine if a domain is a surface.
% a = isSurface(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     1 if the domain is a surface, and 0 otherwise.
%

t = domainType(d);
if t == 102
    a = 1;
else
    a = 0;
end
