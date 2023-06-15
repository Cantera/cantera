function s = importInterface(file, name, phase1, phase2)
% IMPORTINTERFACE  Import an interface between phases.
% s = importInterface(file, name, phase1, phase2)
%
% See `ideal-surface <https://cantera.org/documentation/docs-2.6/sphinx/html/yaml/phases.html#sec-yaml-ideal-surface>`__
% and `Declaring adjacent phases <https://cantera.org/tutorials/yaml/phases.html#declaring-adjacent-phases>`__.
%
% :param file:
%     YAML file containing the interface
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
