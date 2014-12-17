function ok = ready(w)
% READY  Check whether a wall is ready.
% ok = ready(w)
% :param w:
%     Instance of class :mat:func:`Wall`
% :return:
%     Status of the wall
%

ok = wallmethods(11, wall_hndl(w));
