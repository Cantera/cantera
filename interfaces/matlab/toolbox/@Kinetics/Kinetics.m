function k = Kinetics(ph, src, id, neighbor1, neighbor2, neighbor3, neighbor4)
% KINETICS  Kinetics class constructor.
% k = Kinetics(r, ph, neighbor1, neighbor2, neighbor3, neighbor4)
% Class Kinetics represents kinetics managers, which are classes
% that manage reaction mechanisms.  The reaction mechanism
% attributes are specified in a CTML file.
% Instances of class :mat:func:`Kinetics` are responsible for evaluating reaction rates
% of progress, species production rates, and other quantities pertaining to
% a reaction mechanism.
%
% :param ph:
%     An instance of class :mat:func:`ThermoPhase` representing the phase
%     in which reactions occur
% :param src:
%     Input string of YAML, CTI, or XML file name.
% :param id:
%     ID of the phase to import as specified in the input file. (optional)
% :param neighbor1:
%     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
%     neighboring phase.
% :param neighbor2:
%     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
%     neighboring phase.
% :param neighbor3:
%     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
%     neighboring phase.
% :param neighbor4:
%     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
%     neighboring phase.
% :return:
%      Instance of class :mat:func:`Kinetics`
%

% indices for bulk phases in a heterogeneous mechanism.
% initialize < 0 so that bulk phases will not be included.
ineighbor1 = -1;
ineighbor2 = -1;
ineighbor3 = -1;
ineighbor4 = -1;
if nargin == 2
    id = '-';
end

k.owner = 1;
% get the integer indices used to find the stored objects
% representing the phases participating in the mechanism.
iphase = thermo_hndl(ph);
if nargin > 3
    ineighbor1 = thermo_hndl(neighbor1);
    if nargin > 4
        ineighbor2 = thermo_hndl(neighbor2);
        if nargin > 5
            ineighbor3 = thermo_hndl(neighbor3);
            if nargin > 6
                ineighbor4 = thermo_hndl(neighbor4);
            end
        end
    end
end
k.id = kinetics_get(0, 0, src, id, iphase, ineighbor1, ineighbor2, ...
    ineighbor3, ineighbor4);
if k.id < 0
    error(geterr);
end

k = class(k, 'Kinetics');

