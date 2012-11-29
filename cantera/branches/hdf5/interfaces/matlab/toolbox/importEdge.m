function s = importEdge(file, name, phase1, phase2, phase3, phase4)
% IMPORTINTERFACE - import an interface
%
if nargin == 3
    s = Interface(file, name, phase1);
elseif nargin == 4
    s = Interface(file, name, phase1, phase2);
elseif nargin == 5
    s = Interface(file, name, phase1, phase2, phase3);
elseif nargin == 6
    s = Interface(file, name, phase1, phase2, phase3, phase4);
else
    error('importEdge only supports 4 neighbor phases');
end
