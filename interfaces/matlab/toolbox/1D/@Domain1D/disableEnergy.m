function d = disableEnergy(d)
% DISABLEENERGY  Disable the energy equation.
% d = disableEnergy(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
%

disp(' ');
disp('Disabling the energy equation...');

domain_methods(d.dom_id, 66, 0);
