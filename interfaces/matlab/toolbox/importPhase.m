function s = importPhase(file, name)
% IMPORTPHASE  Import a phase from a CTI file
% s = importPhase(file, name)
% Deprecated. To be removed after Cantera 2.3.
% See :ref:`sec-phases`.
%
% See also: :mat:func:`Solution`
%
% :param file:
%     CTI file containing phase definition
% :param name:
%     Name of the phase
% :return:
%     Instance of class :mat:func:`Solution`
%

warning('This function is deprecated and will be removed after Cantera 2.3. Use Solution instead');
if nargin == 1
    s = Solution(file);
elseif nargin == 2
    s = Solution(file, name);
end
