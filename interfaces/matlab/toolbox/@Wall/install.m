function install(w, left, right)
% INSTALL  Install a wall between two reactors.
% install(w, left, right)
% :param w:
%     Instance of class :mat:func:`Wall`
% :param left:
%     Instance of class :mat:func:`Reactor`or
%     :mat:func:`Reservoir`
% :param right:
%     Instance of class :mat:func:`Reactor` or
%     :mat:func:`Reservoir`
%

w.left = left;
w.right = right;
wallmethods(4, wall_hndl(w), reactor_hndl(left), reactor_hndl(right));
