function setMassFractions(tp, y, norm)
% SETMASSFRACTIONS  Set the species mass fractions.
% setMassFractions(tp, y, norm)
% Note that calling :mat:func:`setMassFractions` leaves the temperature and
% density unchanged, and therefore the pressure changes if the new
% composition has a molar mass that is different than the old
% composition. If it is desired to change the composition and hold
% the pressure fixed, use method :mat:func:`set` and specify the mass
% fractions and the pressure, or call :mat:func:`setPressure`
% after calling :mat:func:`setMassFractions`.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param y:
%     Vector of mass fractions whose length must be the same as
%     the number of species. Alternatively, a string in the format
%     ``'SPEC:Y,SPEC2:Y2'`` that specifies the mass fraction of
%     specific species.
% :param norm:
%     If ``'nonorm'`` is specified, ``y`` will be normalized. This only
%     works if ``y`` is a vector, not a string. Since unnormalized mass
%     fractions can lead to unphysical results, ``'nonorm'`` should be
%     used only in rare cases, such as computing partial
%     derivatives with respect to a species mass fraction.
%

if isa(y, 'double')
    if nargin == 3
        if strcmp(norm, 'nonorm')
            phase_set(tp.tp_id, 23, y);
        else
            phase_set(tp.tp_id, 21, y);
        end
    else
        phase_set(tp.tp_id, 21, y);
    end
    %
    % string input
    %
elseif isa(y, 'char')
    phase_set(tp.tp_id, 31, y);
end
