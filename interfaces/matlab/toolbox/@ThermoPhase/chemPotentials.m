function mu = chemPotentials(tp)
% CHEMPOTENTIALS  Get the chemical potentials of the species.
% mu = chemPotentials(tp)
% The expressions used to compute the chemical potential
% depend on the model implemented by the underlying kernel
% thermo manager.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase).
% :return:
%     Vector of species chemical potentials. Units: J/kmol
%
%        This method returns an array containing the species
%        chemical potentials [J/kmol]. The expressions used to
%        compute these depend on the model implemented by the
%        underlying kernel thermo manager."""

mu = thermo_get(tp.tp_id, 34);

