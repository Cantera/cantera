function clear(tr)
% CLEAR  Delete the Transport instance.
% clear(tr)
% :param tr:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%

trans_get(tr.id, 0)
