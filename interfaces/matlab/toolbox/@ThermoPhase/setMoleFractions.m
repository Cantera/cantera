function setMoleFractions(tp, x, norm)
% SETMOLEFRACTIONS  Set the species mole fractions.
% setMoleFractions(tp,x,norm)
% Note that calling :mat:func:`setMoleFractions` leaves the temperature and
% density unchanged, and therefore the pressure changes if the new
% composition has a molar mass that is different than the old
% composition. If it is desired to change the composition and hold
% the pressure fixed, use method :mat:func:`set` and specify the mole
% fractions and the pressure, or call :mat:func:`setPressure`
% after calling :mat:func:`setmoleFractions`.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param y:
%     Vector of mole fractions whose length must be the same as
%     the number of species. Alternatively, a string in the format
%     ``'SPEC:Y,SPEC2:Y2'`` that specifies the mole fraction of
%     specific species.
% :param norm:
%     If ``'nonorm'`` is specified, ``y`` will be normalized. This only
%     works if ``y`` is a vector, not a string. Since unnormalized mole
%     fractions can lead to unphysical results, ``'nonorm'`` should be
%     used only in rare cases, such as computing partial
%     derivatives with respect to a species mole fraction.
%

if isa(x, 'double')
    if nargin == 3
        if strcmp(norm, 'nonorm')
            phase_set(tp.tp_id, 22, x);
        else
            phase_set(tp.tp_id, 20, x);
        end
    else
        phase_set(tp.tp_id, 20, x);
    end
    %
    % string input
    %
elseif isa(x, 'char')
    phase_set(tp.tp_id, 30, x);
end
