function s = importInterface(file, name, phase1, phase2)
    % Import an interface between phases.

    if nargin == 3
        s = Interface(file, name, phase1);
    elseif nargin == 4
        s = Interface(file, name, phase1, phase2);
    else
        error('importInterface only supports 2 bulk phases');
    end
end
