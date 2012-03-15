function install(w, left, right)
% INSTALL -
%
w.left = left;
w.right = right;
wallmethods(4, wall_hndl(w), reactor_hndl(left), ...
    reactor_hndl(right));
