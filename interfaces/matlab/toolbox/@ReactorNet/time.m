function t = time(r)
% TIME  Get the current value of the time.
% t = time(r)
% :param r:
%     Instance of class :mat:func:`ReactorNet`
% :return:
%     Current time in the input ReactorNet. Units: s
%

t = reactornetmethods(22, reactornet_hndl(r));
