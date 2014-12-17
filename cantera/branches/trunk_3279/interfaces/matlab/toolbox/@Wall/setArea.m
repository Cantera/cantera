function setArea(w, a)
% SETAREA  Set the area of a wall.
% setArea(w, a)
% :param w:
%     Instance of class :mat:func:`Wall`
% :param a:
%     Double area of the wall.
%

wallmethods(5, wall_hndl(w), a);
