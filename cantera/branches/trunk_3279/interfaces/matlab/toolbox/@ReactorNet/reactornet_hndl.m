function i = reactornet_hndl(r)
% REACTORNET_HNDL  Get the integer used to access the kernel object.
% i = reactornet_hndl(r)
% :param r:
%     Instance of class :mat:func:`ReactorNet`
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object.
%

i = r.index;
