function advanceCoverages(k, dt)
% ADVANCECOVERAGES - advance the surface coverages forward in time holding the bulk phase concentrations fixed.
%
kinetics_set(k.id, 5, 0, dt);
