function d = enableEnergy(d)
% ENABLEENERGY  Enable the energy equation.
% d = enableEnergy(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
%

disp(' ');
disp('Enabling the energy equation...');

domain_methods(d.dom_id, 66, 1);
