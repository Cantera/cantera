function setMaxJacAge(s, ss_age, ts_age)
% SETMAXJACAGE - Set the number of times the Jacobian will be used
% before it is recomputed.
%
if nargin == 2
    ts_age = ss_age;
end
stack_methods(s.stack_id, 114, ss_age, ts_age);
