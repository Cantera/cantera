function t = atol(r)
% ATOL  Get the absolute error tolerance.
% t = atol(r)
% :param r:
%     Instance of class :mat:func:`ReactorNet`
% :return:
%     Absolute error tolerance.
%

t = reactornetmethods(24, reactornet_hndl(r));

