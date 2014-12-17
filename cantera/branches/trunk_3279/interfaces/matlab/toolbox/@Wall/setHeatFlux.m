function setHeatFlux(w, f)
% SETHEATFLUX  Set the heat flux using :mat:func:`Func`
% setHeatFlux(w, f)
% Must be set by an instance of
% class :mat:func:`Func`, which allows the heat flux to be an
% arbitrary function of time. It is possible to specify
% a constant heat flux by using the polynomial functor with
% only the first term specified:
%
% .. code-block:: matlab
%
%     w = Wall()
%     f = Func('polynomial',0,10); % Or f = polynom(10);
%     setHeatFlux(w, f);
%
% sets the heat flux through the wall to 10 W/m**2.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param f:
%     Instance of class :mat:func:`Func`. Units: W/m**2
%

wallmethods(8, wall_hndl(w), func_hndl(f));
