function a = area(s)
% AREA  Get the area of the reactor surface.
% a = area(s)
% :param s:
%     Instance of class :mat:func:`ReactorSurface`
% :return:
%     Area of the reactor surface in m**2
%

a = reactorsurfacemethods(23, reactorsurface_hndl(s));
