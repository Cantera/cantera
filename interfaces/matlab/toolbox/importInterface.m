function s = importInterface(file, name, phase1, phase2)
% IMPORTINTERFACE  Import an interface between phases.
% s = importInterface(file, name, phase1, phase2)
% See :ref:`sec-interfaces`.
%
% :param file:
%     CTI or CTML file containing the interface
% :param name:
%     Name of the interface to import
% :param phase1:
%     First phase in the interface
% :param phase2:
%     Second phase in the interface
% :return:
%     Instance of class :mat:func:`Interface`
%

if nargin == 3
    s = Interface(file, name, phase1);
elseif nargin == 4
    s = Interface(file, name, phase1, phase2);
else
    error('importInterface only supports 2 bulk phases');
end
