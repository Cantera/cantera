function mu = chemPotentials(p)
% CHEMPOTENTIALS - Species chemical potentials.
%
%        This method returns an array containing the species
%        chemical potentials [J/kmol]. The expressions used to
%        compute these depend on the model implemented by the
%        underlying kernel thermo manager."""
mu = thermo_get(p.tp_id,34);

