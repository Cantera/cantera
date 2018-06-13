function setArea(s, a)
% SETAREA  Set the area of a reactor surface.
% setArea(s, a)
% :param s:
%     Instance of class :mat:func:`ReactorSurface`
% :param a:
%     Double area of the reactor surface.
%

reactorsurfacemethods(5, reactorsurface_hndl(s), a);
