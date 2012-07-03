function x = atomicMasses(a)
% ATOMICMASSES - Array of element atomic masses [kg/kmol].
%
x = phase_get(a.tp_id,30);
