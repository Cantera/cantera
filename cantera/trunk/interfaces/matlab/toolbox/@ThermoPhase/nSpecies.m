function n = nSpecies(tp)
% NSPECIES - Number of species in the phase.
%

n = phase_get(tp.tp_id, 11);
