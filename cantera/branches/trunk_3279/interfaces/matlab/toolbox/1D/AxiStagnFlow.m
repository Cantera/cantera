function m = AxiStagnFlow(gas)
% AXISTAGNFLOW  Get an axisymmetric stagnation flow domain.
% m = AxiStagnFlow(gas)
% :param gas:
%     Instance of class :mat:func:`Solution`
% :return:
%     Domain1D instance representing an axisymmetric
%     stagnation flow.

m = Domain1D(1, gas);
