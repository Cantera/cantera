function restore(s, fname, id)
% RESTORE  Restore a previously-saved solution.
% restore(s, fname, id)
% This method can be used to provide an initial guess for the solution.
%
% See also: :mat:func:`save`
%
% :param s:
%     Instance of class :mat:func:`Stack`
% :param fname:
%     File name of a YAML file containing solution information
% :param id:
%     ID of the element that should be restored
%

stack_methods(s.stack_id, 111, fname, id);
