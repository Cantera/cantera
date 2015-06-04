function v = vdot(w, t)
% VDOT  Get the rate of volumetric change at a given time.
% v = vdot(w, t)
% A positive value corresponds to the left-hand reactor volume
% increasing, and the right-hand reactor volume decreasing.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param t:
%     Time at which the volumetric change should be calculated.
% :return:
%     Rate of volumetric change Units: m**3/s
%

v = wallmethods(21, wall_hndl(w), t);
