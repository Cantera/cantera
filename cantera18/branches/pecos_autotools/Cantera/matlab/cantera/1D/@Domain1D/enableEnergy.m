function d = enableEnergy(d)
% ENABLEENERGY - enable the energy equation
%
disp(' ');
disp('Enabling the energy equation...');

domain_methods(d.dom_id, 66, 1);
