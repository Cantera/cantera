function a = setPressure(tp, p)
% SETPRESSURE  Set the pressure [Pa].
%
%    The pressure  is set by changing the density holding the
%    temperature and chemical composition fixed.
%

if p <= 0.0
    error('The pressure must be positive.')
end

thermo_set(tp.tp_id, 1, p);
