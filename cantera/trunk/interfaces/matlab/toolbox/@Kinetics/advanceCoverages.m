function advanceCoverages(k, dt)
% ADVANCECOVERAGES  Advance the surface coverages forward in time.
% advanceCoverages(k, dt)
% The bulk phase concentrations are held fixed during this operation.
%
% :param k:
%     Instance of class :mat:func:`Interface` with an associated
%     :mat:func:`Kinetics` object.
% :param dt:
%     Time interval by which the coverages should be advanced
%

kinetics_set(k.id, 5, 0, dt);
