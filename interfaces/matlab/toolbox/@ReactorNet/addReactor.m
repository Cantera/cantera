function addReactor(net, reactor)
% ADDREACTOR  Add a reactor to a network.
% addReactor(net, reactor)
% :param net:
%     Instance of class :mat:func:`ReactorNet`
% :param reactor:
%     Instance of class :mat:func:`Reactor`
%

reactornetmethods(4, reactornet_hndl(net), reactor_hndl(reactor));
