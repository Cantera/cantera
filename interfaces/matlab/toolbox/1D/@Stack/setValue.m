function setValue(s, n, comp, localPoint, v)
% SETVALUE  Set the value of a single entry in the solution vector.
% setValue(s, n, comp, localPoint, v)
% Example (assuming ``s`` is an instance of :mat:func:`Stack`)::
%
%     setValue(s, 3, 5, 1, 5.6)
%
% This sets component 5 at the leftmost point (local point 1) in domain 3
% to the value 5.6. Note that the local index always begins at 1
% at the left of each domain, independent of the global index of
% the point, which depends on the location of this domain in the
% stack.
%
% :param s:
%     Instance of class :mat:func:`Stack`
% :param n:
%     Domain number
% :param comp:
%     Component number
% :param localPoint:
%     Local index of the grid point in the domain
% :param v:
%     Value
%

stack_methods(s.stack_id, 100, n, comp, localPoint, v);
