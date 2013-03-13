function gas = air
% AIR - create an object representing air.
%
%    Air is modeled as an ideal gas mixture, and several reactions
%    are defined.  The specification is taked from file air.xml.
%
gas = importPhase('air.cti','air');
