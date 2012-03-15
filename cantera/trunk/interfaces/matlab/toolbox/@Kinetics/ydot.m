function v = ydot(a)
% YDOT - Evaluates wdot_k M_k / (density)
%
v = kinetics_get(a.id,24,0);
