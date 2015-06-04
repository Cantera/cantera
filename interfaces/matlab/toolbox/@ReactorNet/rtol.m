function t = rtol(r)
% RTOL  Get the relative error tolerance.
% t = rtol(r)
% :param r:
%     Instance of class :mat:func:`ReactorNet`
% :return:
%     Relative error tolerance.
%

t = reactornetmethods(23, reactornet_hndl(r));
