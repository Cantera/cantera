function s = importEdge(file, name, phase1, phase2, phase3, phase4)
    % Import edges between phases.
    % s = importEdge(file, name, phase1, phase2, phase3, phase4)
    % Supports up to four neighbor phases. See
    %
    % :param file:
    %     File containing phases
    % :param name:
    %     Name of phase
    % :param phase1:
    %     First neighbor phase
    % :param phase2:
    %     Second neighbor phase
    % :param phase3:
    %     Third neighbor phase
    % :param phase4:
    %     Fourth neighbor phase
    % :return:
    %     Instance of class :mat:func:`Interface`
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
        error('importEdge only supports 4 neighbor phases.');
    end
end
