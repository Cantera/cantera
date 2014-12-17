function i = kinetics_hndl(k)
% KINETICS_HNDL  Get the integer used to access kernel object.
% i = kinetics_hndl(k)
% :param k:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the handle is desired.
% :return:
%     Returns the integer ID of the kinetics kernel object.
%

i = k.id;
