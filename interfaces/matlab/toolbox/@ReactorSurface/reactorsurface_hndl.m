function i = reactorsurface_hndl(s)
% REACTORSURFACE_HNDL  Get the integer used to access the Cantera C++ object.
% i = reactorsurf_hndl(s)
% :param s:
%     Instance of class :mat:func:`ReactorSurface` for which the handle is
%     desired.
% :return:
%     Integer used to access the Cantera C++ object
%

i = s.index;
