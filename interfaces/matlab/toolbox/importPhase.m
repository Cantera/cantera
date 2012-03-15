function s = importPhase(file, name)
% IMPORTPHASE - import a phase from a CTI file
%
if nargin == 1
    s = Solution(file);
elseif nargin == 2
    s = Solution(file, name);
end
