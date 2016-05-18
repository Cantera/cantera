function setKinetics(s, kinetics)
% SETKINETICS  Set the surface reaction mechanisms on a reactor surface.
% setKinetics(w, kinetics)
% :param s:
%     Instance of class :mat:func:`ReactorSurface`
% :param kinetics:
%     Instance of class :mat:func:`Kinetics` (or another object derived from
%     Kinetics) to be used as the kinetic mechanism for this surface. Typically
%     an instance of class :mat:func:`Interface`.
%

ikin = 0;
if isa(kinetics, 'Kinetics')
    ikin = kinetics_hndl(kinetics);
end

reactorsurfacemethods(12, reactorsurface_hndl(s), ikin);
