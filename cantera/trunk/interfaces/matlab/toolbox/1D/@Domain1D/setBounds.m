function setBounds(d, component, lower, upper)
% SETBOUNDS  Set bounds on the solution components.
% d = setBounds(d, component, lower, upper)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param component:
%     String, component to set the bounds on
% :param lower:
%     Lower bound
% :param upper:
%     Upper bound
%
n = componentIndex(d, component);
domain_methods(d.dom_id, 51, n, lower, upper);
