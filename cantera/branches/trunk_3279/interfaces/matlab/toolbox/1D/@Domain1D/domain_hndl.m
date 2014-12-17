function n = domain_hndl(d)
% DOMAIN_HNDL  Get the integer used to access the kernel object.
% n = domain_hndl(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object
%

n = d.dom_id;

