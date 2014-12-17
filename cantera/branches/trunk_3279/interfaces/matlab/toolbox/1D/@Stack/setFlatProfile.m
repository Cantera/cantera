function setFlatProfile(s, domain, comp, v)
% SETFLATPROFILE  Set a component to a value across the entire domain.
% setFlatProfile(s, domain, comp, v)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param domain:
%     Integer ID of the domain
% :param comp:
%     Component to be set
% :param v:
%     Double, value to be set
%

stack_methods(s.stack_id, 102, domain, comp, v);
