function q = qdot(w, t)
% QDOT  Get the total heat transfer through a wall given a time.
% q = qdot(w, t)
% A positive value corresponds to heat flowing from the left-hand
% reactor to the right-hand one.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param t:
%     Time at which the heat transfer should be evaluated.
% :return:
%     Total heat transfer. Units: W
%

q = wallmethods(22, wall_hndl(w), t);
