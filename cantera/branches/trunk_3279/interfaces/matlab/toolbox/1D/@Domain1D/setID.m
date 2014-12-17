function setID(d, id)
% SETID  Set the ID tag for a domain.
% d = setID(d, id)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param id:
%     String ID to assign
%

domain_methods(d.dom_id, 54, id);
