function x = atomicMasses(tp)
% ATOMICMASSES - Array of element atomic masses [kg/kmol].
%

x = phase_get(tp.tp_id, 30);
