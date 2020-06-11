function ch = speciesCharges(tp)
% SPECIESCHARGES  Get the elementary charge of the species.
% x = speciesCharges(tp)
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of species charges. Units: elem. charge
%

ch = phase_get(tp.tp_id, 23);