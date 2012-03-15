function d = disableEnergy(d)
% ENABLEENERGY - enable the energy equation
%
domain_methods(d.dom_id, 66, 0);
