function mmw = meanMolecularWeight(tp)
% MEANMOLECULARWEIGHT  Get the mean molecular weight.
% wtm = meanMolecularWeight(tp)
% The mean molecular weight is the mole-fraction-weighted sum of the
% molar masses of the individual species in the phase.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Scalar double mean molecular weight. Units: kg/kmol
%

mmw = phase_get(tp.tp_id,4);
