function setVelocity(w, f)
% SETVELOCITY  Set the velocity of the wall using :mat:func:`Func`.
% setVelocity(w, f)
% Must be set by an instance of class :mat:func:`Func`, which allows
% the velocity to be an arbitrary function of time. It is possible
% to specify a constant velocity by using the polynomial functor with
% only the first term specified:
%
% .. code-block:: matlab
%
%     w = Wall()
%     f = Func('polynomial',0,10); % Or f = polynom(10);
%     setVelocity(w, f);
%
% sets the velocity of the wall to 10 m/s.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param f:
%     Instance of class :mat:func:`Func`. Units: m/s
%

wallmethods(10, wall_hndl(w), func_hndl(f));
