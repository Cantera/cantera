function setMaxJacAge(s, ss_age, ts_age)
% SETMAXJACAGE  Set the number of times the Jacobian will be used before it is recomputed.
% setMaxJacAge(s, ss_age, ts_age)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param ss_age:
%     Maximum age of the Jacobian for steady state analysis
% :param ts_age:
%     Maximum age of the Jacobian for transient analysis. If
%     not specified, defaults to ``ss_age``.
%

if nargin == 2
    ts_age = ss_age;
end
stack_methods(s.stack_id, 114, ss_age, ts_age);
