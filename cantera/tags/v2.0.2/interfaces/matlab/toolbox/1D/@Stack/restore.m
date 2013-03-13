function restore(s, fname, id)
% RESTORE - Restore a previously-saved solution.
%
%   This method can be used to provide an initial guess for the
%   solution.
stack_methods(s.stack_id, 111, fname, id);
