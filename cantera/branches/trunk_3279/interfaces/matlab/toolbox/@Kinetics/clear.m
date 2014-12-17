function clear(k)
% CLEAR  Delete the Kinetics instance.
% clear(k)
% :param k:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%

kinetics_set(k.id, 3, 0, 0);
